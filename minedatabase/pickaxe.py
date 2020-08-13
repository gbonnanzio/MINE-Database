"""Pickaxe.py: This module generates new compounds from user-specified starting
   compounds using a set of SMARTS-based reaction rules."""
import collections
import csv
import datetime
import hashlib
import itertools
import multiprocessing
import os
import re
import time
from sys import exit
from argparse import ArgumentParser
from copy import deepcopy

from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D

from minedatabase import utils
from minedatabase.databases import MINE
import minedatabase.databases as databases
from minedatabase.utils import rxn2hash, StoichTuple, get_fp


class Pickaxe:
    """This class generates new compounds from user-specified starting
    compounds using a set of SMARTS-based reaction rules. It may be initialized
    with a text file containing the reaction rules and coreactants or this may
    be done on an ad hoc basis."""
    def __init__(self, rule_list=None, coreactant_list=None, explicit_h=True,
                 kekulize=True, neutralise=True, errors=True,
                 racemize=False, database=None, database_overwrite=False,
                 mongo_uri='mongodb://localhost:27017',
                 image_dir=None, quiet=False):
        """
        :param rule_list: Path to a list of reaction rules in TSV form
        :type rule_list: str
        :param coreactant_list: Path to list of coreactants in TSV form
        :type coreactant_list: str
        :param explicit_h: Explicitly represent bound hydrogen atoms
        :type explicit_h: bool
        :param kekulize: Kekulize structures before applying reaction rules
        :type kekulize: bool
        :param neutralise: Remove charges on structure before applying reaction
            rules
        :type neutralise: bool
        :param errors: Print underlying RDKit warnings and halt on error
        :type errors: bool
        :param racemize: Enumerate all possible chiral forms of a molecule if
            unspecified stereocenters exist
        :type racemize: bool
        :param database: Name of desired Mongo Database
        :type database: str
        :param database_overwrite: Force overwrite of existing DB
        :type database_overwrite: bool
        :param mongo_uri: URI of mongo deployment
        :type mongo_uri: str
        :param image_dir: Path to desired image folder
        :type image_dir: str
        :param quiet: Silence unbalenced reaction warnings
        :type quiet: bool
        """
        self.operators = {}
        self.coreactants = {}
        self._raw_compounds = {}
        self.compounds = {}
        self.reactions = {}
        self.generation = 0
        self.explicit_h = explicit_h
        self.kekulize = kekulize
        self.racemize = racemize
        self.neutralise = neutralise
        self.image_dir = image_dir
        self.errors = errors
        self.quiet = quiet
        self.fragmented_mols = False
        self.radical_check = False
        self.structure_field = None        
        # For tanimoto filtering
        self.target_fps = []
        self.target_smiles = []
        self.crit_tani = None
        self.tani_filter = False
        # database info
        self.mongo_uri = mongo_uri
        
        print(f'----------------------------------------')
        print("Intializing pickaxe object") 
        # Determine if the specified database is legal
        if database:
            db = MINE(database, self.mongo_uri)
            if database in db.client.list_database_names():
                if database_overwrite:
                    # If db exists, remove db from all of core compounds and drop db
                    print(f"Database {database} already exists. ",
                            "Deleting database and removing from core compound mines.")
                    db.core_compounds.update_many({}, {'$pull' : {'MINES' : database}})
                    db.client.drop_database(database)
                    self.mine = database
                else:
                    print(f"Warning! Database {database} already exists."
                            "Specify database_overwrite as true to delete old database and write new.")
                    exit("Exiting due to database name collision.")
                    self.mine = None
            else:
                self.mine = database
            del(db)
        else:
            self.mine = None   

        # Use RDLogger to catch errors in log file. SetLevel indicates mode (
        # 0 - debug, 1 - info, 2 - warning, 3 - critical). Default is no errors
        logger = RDLogger.logger()
        if not errors:
            logger.setLevel(4)

        # Load coreactants (if any) into Pickaxe object
        if coreactant_list:
            with open(coreactant_list) as infile:
                for coreactant in infile:
                    self._load_coreactant(coreactant)

        # Load rules (if any) into Pickaxe object
        if rule_list:
            self._load_operators(rule_list)   

        print("\nDone intializing pickaxe object")
        print(f'----------------------------------------\n') 

    def load_target_compounds(self, target_compound_file=None, crit_tani=0, 
                                structure_field=None, id_field='id'):
        """
        Loads the target list into an list of fingerprints to later compare to compounds to determine 
        if those compounds should be expanded.

        :param target_compound_file: Path to a file containing compounds as tsv
        :type target_compound_file: basestring
        :param crit_tani: The critical tanimoto cutoff for expansion
        :type crit_tani: float
        :param structure_field: the name of the column containing the
            structure incarnation as Inchi or SMILES (Default:'structure')
        :type structure_field: str
        :param id_field: the name of the column containing the desired
            compound ID (Default: 'id)
        :type id_field: str
        :return: compound SMILES
        :rtype: list
        """
        # TODO: load targets from mongo db
        # Specify options for tanimoto filtering
        self.tani_filter = True
        self.crit_tani = crit_tani

        # Set structure field to None otherwise value determined by
        # load_structures can interfere
        self.structure_field = None

        # Load target compounds
        if target_compound_file:
            for line in utils.file_to_dict_list(target_compound_file):
                mol = self._mol_from_dict(line, structure_field)
                if not mol:
                    continue
                # Add compound to internal dictionary as a target
                # compound and store SMILES string to be returned
                smi = AllChem.MolToSmiles(mol, True)
                c_id = line[id_field]
                # Do not operate on inorganic compounds
                if 'C' in smi or 'c' in smi:
                    AllChem.SanitizeMol(mol)
                    self._add_compound(c_id, smi, 'Target Compound', mol)
                    self.target_smiles.append(smi)                    
                    # Generate fingerprints for tanimoto filtering
                    # TODO: is there a faster way of doing this?
                    fp = AllChem.RDKFingerprint(mol)
                    self.target_fps.append(fp)
        else:
            raise ValueError("No input file or database specified for "
                             "target compounds")

        print(f"{len(self.target_smiles)} target compounds loaded")

        return self.target_smiles    

    def load_compound_set(self, compound_file=None, id_field='id'): 
            """If a compound file is provided, this function loads the compounds
            into its internal dictionary. 

            :param compound_file: Path to a file containing compounds as tsv
            :type compound_file: basestring            
            :param id_field: the name of the column containing the desired
                compound ID (Default: 'id)
            :type id_field: str

            :return: compound SMILES
            :rtype: list
            """
            # TODO: add in database/mongo_uri
            # That behavior will allow compounds to be added from multiple sources

            # Set structure field to None otherwise value determined by
            # load_targets can interfere
            self.structure_field = None

            # load compounds
            compound_smiles = []
            if compound_file:
                for line in utils.file_to_dict_list(compound_file):
                    mol = self._mol_from_dict(line, self.structure_field)
                    if not mol:
                        continue
                    # Add compound to internal dictionary as a starting
                    # compound and store SMILES string to be returned
                    smi = AllChem.MolToSmiles(mol, True)
                    c_id = line[id_field]
                    # Do not operate on inorganic compounds
                    if 'C' in smi or 'c' in smi:
                        AllChem.SanitizeMol(mol)
                        self._add_compound(c_id, smi,
                                           cpd_type='Starting Compound', mol=mol)
                        compound_smiles.append(smi)

            # TODO: If a MINE database is being used instead, search for compounds
            # annotated as starting compounds and return those as a list of
            # SMILES strings
            # elif database:
            #     if mongo_uri:
            #         db = MINE(database, mongo_uri)
            #     else:
            #         db = MINE(database, self.mongo_uri)
            #     # Check to see if database has any compounds
            #     if db.compounds.find_one({'Type':'Starting Compound'}):
            #         for compound in db.compounds.find({'Type':'Starting Compound'}):
            #             _id = compound['_id']
            #             smi = compound['SMILES']
            #             # Assume unannotated compounds are starting compounds
            #             if 'type' not in compound:
            #                 compound['Type'] = 'Starting Compound'
            #             self._add_compound(_id, smi, cpd_type=compound['Type'])
            #             compound_smiles.append(smi)
            #     else:
            #         raise ValueError('Specified MINE contains no starting compounds.')
            else:
                raise ValueError('No input file specified for '
                                'starting compounds')

            print(f"{len(compound_smiles)} compounds loaded")
            return compound_smiles

    def transform_compound(self, compound_smiles, rules=None):
        """Perform transformations to a compound returning the products and the
        predicted reactions

        :param compound_smiles: The compound on which to operate represented
            as SMILES
        :type compound_smiles: string
        :param rules: The names of the reaction rules to apply. If none,
            all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as tuple of compound and reaction dicts
        :rtype: tuple
        """
        if not rules:
            rules = self.operators.keys()
        
        # Stage compound for transformation
        # Create Mol object from input SMILES string and remove hydrogens
        # (rdkit)
        mol = AllChem.MolFromSmiles(compound_smiles)
        mol = AllChem.RemoveHs(mol)
        if not mol:
            if self.errors:
                raise ValueError(f"Unable to parse: {compound_smiles}")
            else:
                print(f"Unable to parse: {compound_smiles}")
                return
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)

        # Apply reaction rules to prepared compound
        for rule_name in rules:
            # Lookup rule in dictionary from rule name
            rule = self.operators[rule_name]
            # Get RDKit Mol objects for reactants
            # TODO: add functionality for bimolecular+ reactions
            reactant_mols = tuple([mol if x == 'Any'
                                   else self.coreactants[x][0]
                                   for x in rule[1]['Reactants']])
            # Perform chemical reaction on reactants for each rule
            try:
                product_sets = rule[0].RunReactants(reactant_mols)
            # This error should be addressed in a new version of RDKit
            # TODO: Check this claim
            except RuntimeError:
                print("Runtime ERROR!" + rule_name)
                print(compound_smiles)
                continue
            # No enumeration for reactant stereoisomers. _make_half_rxn
            # returns a generator, the first element of which is the reactants.
            reactants, reactant_atoms = self._make_half_rxn_dev(reactant_mols, rule[1]['Reactants'])
            if not reactants:
                continue
            # By sorting the reactant (and later products) we ensure that
            # compound order is fixed.
            # reactants.sort()
            # TODO: check for balanced rxn and remove
            # Only an issue with Joseph's rules
            for product_mols in product_sets:
                try:
                    # for stereo_prods, product_atoms in self._make_half_rxn(
                    #         product_mols, rule[1]['Products'], self.racemize):
                        # Update predicted compounds list
                    try:
                        products, product_atoms = self._make_half_rxn_dev(product_mols, rule[1]['Products'])
                        if not products:
                            continue
                    except:
                        continue

                    if self._is_atom_balanced(product_atoms, reactant_atoms):
                        # products.sort()
                        for _, cpd_dict in products:
                            if cpd_dict['_id'].startswith('C') and cpd_dict['_id'] not in self.compounds:
                                self.compounds[cpd_dict['_id']] = cpd_dict                        
                        # Get reaction text (e.g. A + B <==> C + D)
                        self._add_reaction_dev(reactants, rule_name,
                                                        products)                          
            
                except (ValueError, MemoryError) as e:
                    if not self.quiet:
                        print(e)
                        print("Error Processing Rule: " + rule_name)
                    continue    

        return self.compounds, self.reactions       
    
    def transform_compound_par(self, cpd_id):
        """Perform transformations to a compound returning the products and the
        predicted reactions

        :param compound_smiles: The compound on which to operate represented
            as SMILES
        :type compound_smiles: string
        :param rules: The names of the reaction rules to apply. If none,
            all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as tuple of compound and reaction dicts
        :rtype: tuple
        """
        if not rules:
            rules = self.operators.keys()
        
        # Stage compound for transformation
        # Create Mol object from input SMILES string and remove hydrogens
        # (rdkit)
        mol = AllChem.MolFromSmiles(compound_smiles)
        mol = AllChem.RemoveHs(mol)
        if not mol:
            if self.errors:
                raise ValueError(f"Unable to parse: {compound_smiles}")
            else:
                print(f"Unable to parse: {compound_smiles}")
                return
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)

        gen_cpds = {}
        gen_rxn = {}
        # Apply reaction rules to prepared compound
        for rule_name in rules:
            # Lookup rule in dictionary from rule name
            rule = self.operators[rule_name]
            # Get RDKit Mol objects for reactants
            # TODO: add functionality for bimolecular+ reactions
            reactant_mols = tuple([mol if x == 'Any'
                                   else self.coreactants[x][0]
                                   for x in rule[1]['Reactants']])
            # Perform chemical reaction on reactants for each rule
            try:
                product_sets = rule[0].RunReactants(reactant_mols)
            # This error should be addressed in a new version of RDKit
            # TODO: Check this claim
            except RuntimeError:
                print("Runtime ERROR!" + rule_name)
                print(compound_smiles)
                continue
            # No enumeration for reactant stereoisomers. _make_half_rxn
            # returns a generator, the first element of which is the reactants.
            reactants, reactant_atoms = self._make_half_rxn_dev(reactant_mols, rule[1]['Reactants'])
            if not reactants:
                continue
            # By sorting the reactant (and later products) we ensure that
            # compound order is fixed.
            # reactants.sort()
            # TODO: check for balanced rxn and remove
            # Only an issue with Joseph's rules
            for product_mols in product_sets:
                try:
                    # for stereo_prods, product_atoms in self._make_half_rxn(
                    #         product_mols, rule[1]['Products'], self.racemize):
                        # Update predicted compounds list
                    try:
                        products, product_atoms = self._make_half_rxn_dev(product_mols, rule[1]['Products'])
                        if not products:
                            continue
                    except:
                        continue

                    if self._is_atom_balanced(product_atoms, reactant_atoms):
                        # products.sort()
                        for _, cpd_dict in products:
                            if cpd_dict['_id'].startswith('C') and cpd_dict['_id'] not in self.compounds:
                                gen_cpds[cpd_dict['_id']] = cpd_dict                        
                        # Get reaction text (e.g. A + B <==> C + D)
                        self._add_reaction_dev(reactants, rule_name,
                                                        products)                          
            
                except (ValueError, MemoryError) as e:
                    if not self.quiet:
                        print(e)
                        print("Error Processing Rule: " + rule_name)
                    continue    

        return self.compounds, self.reactions

    def transform_chunk(self, cpd_id_chunk):
        """This function takes in a chunyuk of compound ids and applies transformations to them,
        storing the new compounds and reactions in a local dictionary to be passed back to,
        the main transform_all function.
        """
        gen_cpds = {}
        gen_rxns = {}

        for cpd_id in cpd_id_chunk:
            cpd_gen_cpds, cpd_gen_rxns = self.transform_compound_par(cpd_id)
            """new_cpds and new_rxns are going to be returned as
            {cpd_id: cpd_dict}
            {rxn_id: rxn_dict}
            Product of will be added later
            """
            gen_cpds.update(cpd_gen_cpds)
            gen_rxns.update(cpd_gen_rxns)
        
        return gen_cpds, gen_rxns

    def transform_all(self, num_workers=1, max_generations=1):
        """This function applies all of the reaction rules to all the compounds
        until the generation cap is reached.

        :param num_workers: The number of CPUs to for the expansion process.
        :type num_workers: int
        :param max_generations: The maximum number of times an reaction rule
            may be applied
        :type max_generations: int
        """
        def print_progress(done, total):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not done % print_on:
                print(f"Generation {self.generation}: {round(done / total * 100)} percent complete")

        while self.generation < max_generations:
            if self.tani_filter == True:
                if not self.target_fps:
                    print(f'No targets to filter for. Terminating expansion.')
                    return None
                # Flag compounds to be expanded
                if type(self.crit_tani) == list:
                    crit_tani = self.crit_tani[self.generation]
                else:
                    crit_tani = self.crit_tani
                print(f"Filtering out tanimoto < {crit_tani}")
                self._filter_by_tani(num_workers=num_workers)
                n_filtered = 0
                n_total = 0
                # TODO better way to record this instead of looping again
                for cpd_dict in self.compounds.values():
                    if cpd_dict['Generation'] == self.generation:
                        n_total += 1
                        if cpd_dict['Expand'] == True:
                            n_filtered += 1

                print(f'{n_filtered} of {n_total} compounds remain after filter.\n\nExpanding.')

            self.generation += 1
            # Use to print out time per generation at end of loop
            time_init = time.time()
            n_comps = len(self.compounds)
            n_rxns = len(self.reactions)
            # Get all SMILES strings for compounds
            compound_smiles = [cpd['SMILES'] for cpd in self.compounds.values()
                            if cpd['Generation'] == self.generation - 1
                            and cpd['Type'] not in ['Coreactant', 'Target Compound']
                            and cpd['Expand'] == True]

            if not compound_smiles:
                continue
            if num_workers > 1:
                chunk_size = max(
                    [round(len(compound_smiles) / (num_workers * 10)), 1])
                print(f"Chunk Size for generation {self.generation}:", chunk_size)
                new_comps = deepcopy(self.compounds)
                new_rxns = deepcopy(self.reactions)
                pool = multiprocessing.Pool(processes=num_workers)
                for i, res in enumerate(pool.imap_unordered(
                        self.transform_compound, compound_smiles, chunk_size)):
                    new_comps.update(res[0])
                    new_rxns.update(res[1])
                    print_progress((i + 1), len(compound_smiles))
                self.compounds = new_comps
                self.reactions = new_rxns
            else:
                for i, smi in enumerate(compound_smiles):
                    # Perform possible reactions on compound
                    try:
                        self.transform_compound(smi)
                        print_progress(i, len(compound_smiles))
                    except:
                        # TODO what error
                        continue

            print(f"Generation {self.generation} took {time.time()-time_init} sec and produced:")
            print(f"\t\t{len(self.compounds) - n_comps} new compounds")
            print(f"\t\t{len(self.reactions) - n_rxns} new reactions")
            print(f'----------------------------------------\n')

    def _load_coreactant(self, coreactant_text):
        """
        Loads a coreactant into the coreactant dictionary from a tab-delimited
            string
        :param coreactant_text: tab-delimited string with the compound name and
            SMILES
        """
        # If coreactant is commented out (with '#') then don't import
        if coreactant_text[0] == '#':
            return
        split_text = coreactant_text.strip().split('\t')
        # split_text[0] is compound name, split_text[1] is SMILES string
        # Generate a Mol object from the SMILES string if possible
        try:
            mol = AllChem.MolFromSmiles(split_text[2])
            if not mol:
                raise ValueError
            # TODO: what do do about stereochemistry? Original comment is below
            # but stereochem was taken out (isn't it removed later anyway?)
            # # Generate SMILES string with stereochemistry taken into account
            smi = AllChem.MolToSmiles(mol)
        except (IndexError, ValueError):
            raise ValueError(f"Unable to load coreactant: {coreactant_text}")
        c_id = self._add_compound(split_text[0], smi, 'Coreactant', mol)
        # If hydrogens are to be explicitly represented, add them to the Mol
        # object
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        # If kekulization is preferred (no aromatic bonds, just 3 C=C bonds
        # in a 6-membered aromatic ring for example)
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        # Store coreactant in a coreactants dictionary with the Mol object
        # and hashed id as values (coreactant name as key)
        self.coreactants[split_text[0]] = (mol, c_id,)

    def _load_operators(self, rule_path):
        """Loads all reaction rules from file_path into rxn_rule dict.

        :param rule_path: path to file
        :type rule_path: str
        """
        skipped = 0
        with open(rule_path) as infile:
            # Get all reaction rules from tsv file and store in dict (rdr)
            rdr = csv.DictReader((row for row in infile if not
                                  row.startswith('#')), delimiter='\t')
            for rule in rdr:
                try:
                    # Get reactants and products for each reaction into list
                    # form (not ; delimited string)
                    rule['Reactants'] = rule['Reactants'].split(';')
                    rule['Products'] = rule['Products'].split(';')
                    # Ensure that all coreactants are known and accounted for
                    all_rules = rule['Reactants'] + rule['Products']
                    for coreactant_name in all_rules:
                        if ((coreactant_name not in self.coreactants
                             and coreactant_name != 'Any')):
                            raise ValueError(f"Undefined coreactant:{coreactant_name}")
                    # Create ChemicalReaction object from SMARTS string
                    rxn = AllChem.ReactionFromSmarts(rule['SMARTS'])
                    rule.update({'_id': rule['Name'],
                                 'Reactions_predicted': 0,
                                 'SMARTS': rule['SMARTS']})
                    # Ensure that we have number of expected reactants for
                    # each rule
                    if rxn.GetNumReactantTemplates() != len(rule['Reactants'])\
                            or rxn.GetNumProductTemplates() != \
                            len(rule['Products']):
                        skipped += 1
                        print("The number of coreactants does not match the "
                              "number of compounds in the SMARTS for reaction "
                              "rule: " + rule['Name'])
                    if rule['Name'] in self.operators:
                        raise ValueError("Duplicate reaction rule name")
                    # Update reaction rules dictionary
                    self.operators[rule['Name']] = (rxn, rule)
                except Exception as e:
                    raise ValueError(str(e) + f"\nFailed to parse {rule['Name']}")
        if skipped:
            print("WARNING: {skipped} rules skipped")  

    def _filter_by_tani(self, num_workers=1):
        """ 
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target compound
        greater than or equal to the crit_tani, for expansion.
        """
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = [cpd for cpd in self.compounds.values() 
                                if cpd['Generation'] == self.generation
                                and cpd['Type'] not in ['Coreactant', 'Target Compound']]
        
        if num_workers > 1:
            # Set up parallel computing of compounds to expand
            chunk_size = max(
                        [round(len(compounds_to_check) / (num_workers * 10)), 1])
            print(f'Filtering Generation {self.generation}')
            pool = multiprocessing.Pool(num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    self._compare_to_targets, 
                    [cpd for cpd in compounds_to_check 
                        if cpd['Generation'] == self.generation], chunk_size)):                
                # If the result of comparison is false, compound is not expanded
                # Default value for a compound is True, so no need to specify expansion
                # TODO: delete these compounds instead of just labeling as false?
                if not res[1]:
                    self.compounds[res[0]]['Expand'] = False
                print_progress(i, len(compounds_to_check), 'Tanimoto filter progress:')
        else:
            print(f'Filtering Generation {self.generation}')
            cpd_to_compare = [cpd for cpd in compounds_to_check 
                        if cpd['Generation'] == self.generation]
            for i, cpd in enumerate(cpd_to_compare):
                res = self._compare_to_targets(cpd)
                if not res[1]:
                    self.compounds[res[0]]['Expand'] = False
                print_progress(i, len(compounds_to_check), 'Tanimoto filter progress:')            
            
        return None
    
    def _compare_to_targets(self, cpd):
        """ 
        Helper function to allow parallel computation of tanimoto filtering.
        Works with _filter_by_tani

        Returns True if a the compound is similar enough to a target.
        
        """
        # Generate the fingerprint of a compound and compare to the fingerprints of the targets
        if type(self.crit_tani) == list:
            crit_tani = self.crit_tani[self.generation]
        else:
            crit_tani = self.crit_tani

        try:
            fp1 = get_fp(cpd['SMILES'])
            for fp2 in self.target_fps:
                if AllChem.DataStructs.FingerprintSimilarity(fp1, fp2) >= crit_tani:
                    return (cpd['_id'], True)
        except:
            pass

        return (cpd['_id'], False)

    def _mol_from_dict(self, input_dict, structure_field=None):
        # detect structure field as needed
        if not structure_field:
            if not self.structure_field:
                for field in input_dict:
                    if str(field).lower() in {'smiles', 'inchi', 'structure'}:
                        self.structure_field = field
                        break
            if not self.structure_field:
                raise ValueError('Structure field not found in input')
            structure_field = self.structure_field

        if structure_field not in input_dict:
            return
        # Generate Mol object from InChI code if present
        if 'InChI=' in input_dict[structure_field]:
            mol = AllChem.MolFromInchi(input_dict[structure_field])
        # Otherwise generate Mol object from SMILES string
        else:
            mol = AllChem.MolFromSmiles(input_dict[structure_field])
        if not mol:
            if self.errors:
                print(f"Unable to Parse {input_dict[structure_field]}")
            return
        # If compound is disconnected (determined by GetMolFrags
        # from rdkit) and loading of these molecules is not
        # allowed (i.e. fragmented_mols == 1), then don't add to
        # internal dictionary. This is most common when compounds
        # are salts.
        if not self.fragmented_mols and len(AllChem.GetMolFrags(mol)) > 1:
            return
        # If specified remove charges (before applying reaction
        # rules later on)
        if self.neutralise:
            mol = utils.neutralise_charges(mol)
        return mol

    def _gen_compound(self, cpd_id, smi, cpd_type, mol=None):
        """Generates a compound"""
        cpd_dict = {}
        c_id = utils.compound_hash(smi, cpd_type)
        self._raw_compounds[smi] = c_id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if c_id not in self.compounds:
            if not mol:
                mol = AllChem.MolFromSmiles(smi)
            # expand only Predicted and Starting_compounds
            expand = True if cpd_type in ['Predicted', 'Starting Compound'] else False
            cpd_dict = {'ID': cpd_id, '_id': c_id, 'SMILES': smi,
                                   'Type': cpd_type,
                                   'Generation': self.generation,
                                   'atom_count': utils._getatom_count(mol, self.radical_check),
                                   'Reactant_in': [], 'Product_of': [],
                                   'Expand': expand}
            if c_id[0] =='X':
                del(cpd_dict['Reactant_in'])
                del(cpd_dict['Product_of'])
        else:
            cpd_dict = self.compounds[c_id]
        
        return c_id, cpd_dict
    
    def _insert_compound(self, cpd_dict):
        """Inserts a compound into the dictionary"""        
        self.compounds[cpd_dict['_id']] = cpd_dict
        # if self.image_dir and self.mine:
        #         try:
        #             with open(os.path.join(self.image_dir, c_id + '.svg'),
        #                     'w') as outfile:
        #                 nmol = rdMolDraw2D.PrepareMolForDrawing(mol)
        #                 d2d = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
        #                 d2d.DrawMolecule(nmol)
        #                 d2d.FinishDrawing()
        #                 outfile.write(d2d.GetDrawingText())
        #         except OSError:
        #             print(f"Unable to generate image for {smi}")

    def _add_compound(self, cpd_id, smi, cpd_type, mol=None):
        """Adds a compound to the internal compound dictionary"""
        c_id = utils.compound_hash(smi, cpd_type)
        self._raw_compounds[smi] = c_id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if c_id not in self.compounds:
            _, cpd_dict = self._gen_compound(cpd_id, smi, cpd_type, mol)
            self._insert_compound(cpd_dict)
        
        return c_id

    def _add_reaction(self, reactants, rule_name, stereo_prods):
        """Hashes and inserts reaction into reaction dictionary"""
        # Hash reaction text
        rhash = rxn2hash(reactants, stereo_prods)
        # Generate unique hash from InChI keys of reactants and products
        # inchi_rxn_hash, text_rxn = \
        #     self._calculate_rxn_hash_and_text(reactants, stereo_prods)
        text_rxn = self._calculate_rxn_text(reactants, stereo_prods)
        # Add reaction to reactions dictionary if not already there
        if rhash not in self.reactions:
            self.reactions[rhash] = {'_id': rhash,
                                     'Reactants': reactants,
                                     'Products': stereo_prods,
                                    #  'InChI_hash': inchi_rxn_hash,
                                     'Operators': {rule_name},                                    
                                     'SMILES_rxn': text_rxn,
                                     'Generation': self.generation}

        # Otherwise, update the operators
        else:
            self.reactions[rhash]['Operators'].add(rule_name)

        # Update compound tracking
        for prod_id in [x.c_id for x in stereo_prods if x.c_id[0] == 'C']:
            if rhash not in self.compounds[prod_id]['Product_of']:
                self.compounds[prod_id]['Product_of'].append(rhash)
        
        for reac_id in [x.c_id for x in reactants if x.c_id[0] == 'C']:
            if rhash not in self.compounds[reac_id]['Reactant_in']:
                self.compounds[reac_id]['Reactant_in'].append(rhash)      

        return text_rxn

    def _add_reaction_dev(self, reactants, rule_name, products):
        """Hashes and inserts reaction into reaction dictionary"""
        # Hash reaction text
        rhash = rxn2hash(reactants, products)
        # Generate unique hash from InChI keys of reactants and products
        # inchi_rxn_hash, text_rxn = \
        #     self._calculate_rxn_hash_and_text(reactants, stereo_prods)
        # STOPPING HERE
        text_rxn = self._calculate_rxn_text(reactants, products)
        # Add reaction to reactions dictionary if not already there
        if rhash not in self.reactions:
            self.reactions[rhash] = {'_id': rhash,
                                     'Reactants': reactants,
                                     'Products': products,
                                    #  'InChI_hash': inchi_rxn_hash,
                                     'Operators': {rule_name},                                    
                                     'SMILES_rxn': text_rxn}

        # Otherwise, update the operators
        else:
            self.reactions[rhash]['Operators'].add(rule_name)

        # Update compound tracking
        for prod_id in [cpd_dict['_id'] for _, cpd_dict in products if cpd_dict['_id'].startswith('C')]:
            if rhash not in self.compounds[prod_id]['Product_of']:
                self.compounds[prod_id]['Product_of'].append(rhash)
        
        for reac_id in [cpd_dict['_id'] for _, cpd_dict in reactants if cpd_dict['_id'].startswith('C')]:
            if rhash not in self.compounds[reac_id]['Reactant_in']:
                self.compounds[reac_id]['Reactant_in'].append(rhash)      

        return text_rxn

    def _add_reaction_par(self, reactants, rule_name, products):
        """Hashes and returns rxn dicts"""
        # Hash reaction text
        rhash = rxn2hash(reactants, products)
        text_rxn = self._calculate_rxn_text(reactants, products)
        # Add reaction to reactions dictionary if not already there
        if rhash not in self.reactions:
            self.reactions[rhash] = {'_id': rhash,
                                     'Reactants': reactants,
                                     'Products': products,
                                    #  'InChI_hash': inchi_rxn_hash,
                                     'Operators': {rule_name},                                    
                                     'SMILES_rxn': text_rxn}

        # Otherwise, update the operators
        else:
            self.reactions[rhash]['Operators'].add(rule_name)

        # Update compound tracking
        for prod_id in [cpd_dict['_id'] for _, cpd_dict in products if cpd_dict['_id'].startswith('C')]:
            if rhash not in self.compounds[prod_id]['Product_of']:
                self.compounds[prod_id]['Product_of'].append(rhash)
        
        for reac_id in [cpd_dict['_id'] for _, cpd_dict in reactants if cpd_dict['_id'].startswith('C')]:
            if rhash not in self.compounds[reac_id]['Reactant_in']:
                self.compounds[reac_id]['Reactant_in'].append(rhash)      

        return text_rxn

    def _is_atom_balanced(self, product_atoms, reactant_atoms):
        """If the SMARTS rule is not atom balanced, this check detects the
        accidental alchemy."""
        if reactant_atoms - product_atoms \
                or product_atoms - reactant_atoms:
            return False
        else:
            return True      

    def _make_half_rxn_dev(self, mol_list, rules):
        """Takes a list of mol objects for a half reaction, combines like
        compounds and returns a generator for stoich tuples"""
        # Get compound ids from Mol objects, except for coreactants, in which
        #  case we look them up in the coreactant dictionary
        cpds = {}
        cpd_counter = collections.Counter()

        for mol, rule in zip(mol_list, rules):
            if rule == 'Any':
                cpd_id, cpd_dict = self._calculate_compound_information_dev(mol)
                if cpd_id == None:
                    # failed to make rxn
                    return None, None
            else:
                cpd_id = self.coreactants[rule][1]
                cpd_dict = self.compounds[cpd_id]

            cpds[cpd_id] = cpd_dict
            cpd_counter.update({cpd_id : 1})

        # count the number of atoms on a side
        atom_counts = collections.Counter()
        for cpd_id, cpd_dict in cpds.items():
            # atom_count += cpd_dict['atom_count']*cpd_counter[cpd_id]        
            for atom_id, atom_count in cpd_dict['atom_count'].items():
                atom_counts[atom_id] += atom_count*cpd_counter[cpd_id]

        return [(stoich, cpds[cpd_id]) for cpd_id, stoich in cpd_counter.items()], atom_counts

    def _calculate_compound_information_dev(self, mol_obj):
        """Calculate the standard data for a compound & return a tuple with
        compound_ids. Memoized with _raw_compound dict"""
        # This is the raw SMILES which may have explicit hydrogen
        raw = AllChem.MolToSmiles(mol_obj, True)
        if raw not in self._raw_compounds:
            try:
                # Remove hydrogens if explicit (this step slows down the
                # process quite a bit)
                if self.explicit_h:
                    mol_obj = AllChem.RemoveHs(mol_obj)
                AllChem.SanitizeMol(mol_obj)
            except:
                return None, None
            # In case we want to have separate entries for stereoisomers
            
            # Get the molecule smiles
            smiles = AllChem.MolToSmiles(mol_obj, True)
            cpd_id, cpd_dict = self._gen_compound(None, smiles, 'Predicted', mol_obj)
            self._raw_compounds[raw] = cpd_id
        else:
            cpd_id = self._raw_compounds[raw]
            cpd_dict = self.compounds[cpd_id]
            

        return cpd_id, cpd_dict

    def _calculate_rxn_text(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""
        def get_blocks(cpds):
            cpd_tups = [(stoich, cpd_dict['_id'], cpd_dict['SMILES']) for stoich, cpd_dict in cpds]
            cpd_tups.sort(key=lambda x: x[1])
            smiles = []
            for cpd in cpd_tups:
                smiles.append(f"({cpd[0]}) {cpd[2]}")                
            return ' + '.join(smiles)

        r_s = get_blocks(reactants)
        p_s = get_blocks(products)
        smiles_rxn = r_s + ' => ' + p_s        
        return smiles_rxn

    def _calculate_rxn_hash_and_text(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""
        def __get_blocks(tups):
            first_block, second_block, smiles = [], [], []
            for x in tups:
                comp = self.compounds[x.c_id]
                smiles.append(f"({x.stoich}) {comp['SMILES']}")
                if comp['Inchikey']:
                    # InChI keys are separated by a hyphen, where the first
                    # part is derived from connectivity and the second part
                    # comes from other layers such as stereochemistry
                    split_inchikey = comp['Inchikey'].split('-')
                    if len(split_inchikey) > 1:
                        first_block.append(f'{x.stoich},{split_inchikey[0]}')
                        second_block.append(f'{x.stoich},{split_inchikey[1]}')
                else:
                    print(f"No Inchikey for {x.c_id}")
            return '+'.join(first_block), '+'.join(second_block), \
                   ' + '.join(smiles)

        reactants.sort()
        products.sort()
        r_1, r_2, r_s = __get_blocks(reactants)
        p_1, p_2, p_s = __get_blocks(products)
        first_block = r_1 + '<==>' + p_1
        second_block = r_2 + '<==>' + p_2
        smiles_rxn = r_s + ' => ' + p_s

        return hashlib.sha256(first_block.encode()).hexdigest() + '-' \
            + hashlib.md5(second_block.encode()).hexdigest(), smiles_rxn

    def assign_ids(self):
        """Assigns a numerical ID to compounds (and reactions) for ease of
        reference. Unique only to the CURRENT run."""
        # If we were running a multiprocess expansion, this removes the dicts
        # from Manager control
        self.compounds = dict(self.compounds)
        self.reactions = dict(self.reactions)
        i = 1
        for comp in sorted(self.compounds.values(),
                           key=lambda x: (x['Generation'], x['_id'])):
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            if not comp['ID']:
                comp['ID'] = 'pkc' + str(i).zfill(7)
                i += 1
                self.compounds[comp['_id']] = comp
                # If we are not loading into the mine, we generate the image
                # here.
                if self.image_dir and not self.mine:
                    mol = AllChem.MolFromSmiles(comp['SMILES'])
                    try:
                        MolToFile(
                            mol,
                            os.path.join(self.image_dir, comp['ID'] + '.png'),
                            fitImage=True, kekulize=False)
                    except OSError:
                        print(f"Unable to generate image for {comp['SMILES']}")
        i = 1
        for rxn in sorted(self.reactions.values(),
                          key=lambda x: (x['Generation'], x['_id'])):
            rxn['ID_rxn'] = ' + '.join(
                [f"({x.stoich}) {self.compounds[x.c_id]['ID']}[c0]"
                 for x in rxn['Reactants']]) + ' => ' + ' + '.join(
                     [f"({x.stoich}) {self.compounds[x.c_id]['ID']}[c0]"
                      for x in rxn['Products']])
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            rxn['ID'] = 'pkr' + str(i).zfill(7)
            i += 1
            self.reactions[rxn['_id']] = rxn

    def prune_network_to_targets(self):
        """
        Remove compounds that were unexpanded as well as reactions that ended terminally with them.
        """
        white_list = []
        for target_smi in self.target_smiles:
            try:
                # generate hash of predicted target compounds
                white_list.append(utils.compound_hash(target_smi, 'Predicted'))
            except:
                pass
        
        self.prune_network(white_list)

    def prune_network(self, white_list):
        """
        Prune the predicted reaction network to only compounds and reactions
        that terminate in a specified white list of compounds.
        :param white_list: A list of compound_ids to include (if found)
        :type white_list: list
        :return: None
        """
        n_white = len(white_list)
        comp_set, rxn_set = self.find_minimal_set(white_list)
        print(f"Pruned network to {len(comp_set)} compounds and {len(rxn_set)} reactions based on \
                {n_white} whitelisted compounds")
        self.compounds = dict([(k, v) for k, v in self.compounds.items()
                               if k in comp_set])
        self.reactions = dict([(k, v) for k, v in self.reactions.items()
                               if k in rxn_set])

    def find_minimal_set(self, white_list):
        """
        Given a whitelist this function finds the minimal set of compound and
        reactions ids that comprise the set
        :param white_list:  A list of compound_ids to include (if found)
        :type white_list: list
        :return: compound and reaction id sets
        :rtype: tuple(set, set)
        """
        white_set = set(white_list)
        comp_set = set()
        rxn_set = set()
        for c_id in white_list:
            if c_id not in self.compounds:
                continue
            for r_id in self.compounds[c_id]['Product_of']:
                rxn_set.add(r_id)
                comp_set.update([x.c_id for x
                                 in self.reactions[r_id]['Products']])
                for tup in self.reactions[r_id]['Reactants']:
                    comp_set.add(tup.c_id)
                    # do not want duplicates or cofactors in the whitelist
                    if tup.c_id[0] == 'C' and tup.c_id not in white_set:
                        white_list.append(tup.c_id)
                        white_set.add(tup.c_id)

        # Save targets
        # TODO: seperate targets
        if self.tani_filter:
           for c_id in self.compounds:
               if c_id.startswith('T'):
                   comp_set.add(c_id)

        return comp_set, rxn_set

    def write_compound_output_file(self, path, dialect='excel-tab'):
        """Writes all compound data to the specified path.

        :param path: path to output
        :type path: str
        :param dialect: the output format for the file. Choose excel for csv
            excel-tab for tsv.
        :type dialect: str
        """
        path = utils.prevent_overwrite(path)

        columns = ('ID', 'Type', 'Generation', 'Formula', 'Inchikey',
                   'SMILES')
        with open(path, 'w') as outfile:
            writer = csv.DictWriter(outfile, columns, dialect=dialect,
                                    extrasaction='ignore', lineterminator='\n')
            writer.writeheader()
            writer.writerows(sorted(self.compounds.values(),
                                    key=lambda x: x['ID']))

    def write_reaction_output_file(self, path, delimiter='\t'):
        """Writes all reaction data to the specified path.

        :param path: path to output
        :type path: basestring
        :param delimiter: the character with which to separate data entries
        :type delimiter: basestring
        """
        path = utils.prevent_overwrite(path)
        with open(path, 'w') as outfile:
            outfile.write('ID\tName\tID equation\tSMILES equation\tRxn hash\t'
                          'Reaction rules\n')
            for rxn in sorted(self.reactions.values(), key=lambda x: x['ID']):
                outfile.write(delimiter.join([rxn['ID'], '', rxn['ID_rxn'],
                                              rxn['SMILES_rxn'], rxn['_id'],
                                              ';'.join(rxn['Operators'])])
                              + '\n')
    
    def _save_compounds(self, cpd_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then delete them from memory
        core_cpd_requests = []
        core_update_mine_requests = []
        mine_cpd_requests = []

        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(cpd_ids) / (num_workers * 10)), 1])
            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    self._save_compound_helper, cpd_ids, chunk_size)):
                if res:
                    mine_cpd_requests.append(res[0])
                    core_update_mine_requests.append(res[1])
                    core_cpd_requests.append(res[2])    
        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for cpd_id in cpd_ids:
                comp_dict = self.compounds[cpd_id]
                # Write everything except for targets
                if not comp_dict['_id'].startswith('T'):
                    # These functions are in the MINE class. The request list is
                    # passed and appended in the MINE method.
                    db.insert_mine_compound(comp_dict, mine_cpd_requests)
                    db.update_core_compound_MINES(comp_dict, core_update_mine_requests)
                    db.insert_core_compound(comp_dict, core_cpd_requests)
        return core_cpd_requests, core_update_mine_requests, mine_cpd_requests
    
    def _save_compound_helper(self, cpd_id):
        # Helper function to aid parallelization of saving compounds in
        # save_to_mine
        comp_dict = self.compounds[cpd_id]
        if not comp_dict['_id'].startswith('T'):
            # These functions are outside of the MINE class in order to
            # allow for parallelization. When in the MINE class it is not
            # serializable with pickle. In comparison to the class functions,
            # these return the requests instead of appending to a passed list.
            mine_req = databases.insert_mine_compound(comp_dict)
            core_up_req = databases.update_core_compound_MINES(comp_dict, self.mine)
            core_in_req = databases.insert_core_compound(comp_dict)
            return [mine_req, core_up_req, core_in_req]
        else:
            return None
    
    def _save_reactions(self, rxn_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then delete them from memory
        mine_rxn_requests = []
        rxns_to_write = [self.reactions[rxn_id] for rxn_id in rxn_ids]

        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(rxn_ids) / (num_workers * 10)), 1])
            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    databases.insert_reaction, rxns_to_write, chunk_size)):
                if res:
                    mine_rxn_requests.append(res)

        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for rxn_id in rxn_ids:
                rxn = self.reactions[rxn_id]
                db.insert_reaction(rxn, requests=mine_rxn_requests)
        
        return mine_rxn_requests

    def _save_target_helper(self, comp_dict):
        # Helper function to aid parallelization of saving targets in
        # save_to_mine
        if comp_dict['_id'].startswith('T'):
            # This functions are outside of the MINE class in order to
            # allow for parallelization. When in the MINE class it is not
            # serializable with pickle. In comparison to the class functions,
            # these return the requests instead of appending to a passed list.
            target_req = databases.insert_mine_compound(comp_dict)
            return target_req
        else:
            return None
    
    def save_to_mine(self, num_workers=1, indexing=True):
        """Save compounds to a MINE database.

        :param db_id: The name of the target database
        :type db_id: basestring
        """        
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            n = max(n, 1)           
            for i in range(0, len(lst), n):
                yield lst[i:i + n]
            
        print(f'----------------------------------------')
        print(f'Saving results to {self.mine}')
        print(f'----------------------------------------\n')
        start = time.time()
        db = MINE(self.mine, self.mongo_uri)    

        # Insert Reactions    
        print('--------------- Reactions ---------------')
        rxn_start = time.time()
        # Due to memory concerns, compounds are chunked
        # and processed that way. Each batch is calculated
        # in parallel.
        n_rxns = len(self.reactions)
        chunk_size = max(int(n_rxns/100), 10000)     
        print(f"Reaction chunk size writing: {chunk_size}")   
        n_loop = 1
        for rxn_id_chunk in chunks(list(self.reactions.keys()), chunk_size):
            print(f"Writing Reaction Chunk {n_loop} of {round(n_rxns/chunk_size)+1}")
            n_loop += 1
            mine_rxn_requests = self._save_reactions(rxn_id_chunk, db, num_workers)
            if mine_rxn_requests:
                db.reactions.bulk_write(mine_rxn_requests, ordered=False)
        print(f'Finished Inserting Reactions in {time.time() - rxn_start} seconds.')
        print(f'----------------------------------------\n')

        print('--------------- Compounds --------------')
        cpd_start = time.time()
        # Due to memory concerns, compounds are chunked
        # and processed that way. Each batch is calculated
        # in parallel.
        n_cpds = len(self.compounds)
        chunk_size = max(int(n_cpds/100), 1000)
        print(f"Compound chunk size writing: {chunk_size}")  
        n_loop = 1
        for cpd_id_chunk in chunks(list(self.compounds.keys()), chunk_size):
            # Insert the three types of compounds
            print(f"Writing Compound Chunk {n_loop} of {round(n_cpds/chunk_size)+1}")
            n_loop += 1
            core_cpd_requests, core_update_mine_requests, mine_cpd_requests = self._save_compounds(cpd_id_chunk, db, num_workers)            
            if core_cpd_requests:
                db.core_compounds.bulk_write(core_cpd_requests, ordered=False)            
                del(core_cpd_requests)            
            if core_update_mine_requests:
                db.core_compounds.bulk_write(core_update_mine_requests, ordered=False)            
                del(core_update_mine_requests)
            if mine_cpd_requests:
                db.compounds.bulk_write(mine_cpd_requests, ordered=False)                     
                del(mine_cpd_requests) 
        print(f" Finished insertin Compounds to the MINE in {time.time() - cpd_start} seconds.")
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                        "Action": "Core Compounds Inserted"})
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
            "Action": "Mine Compounds Inserted"})
        print(f"----------------------------------------\n")

        # Insert target compounds
        
        if self.tani_filter:
            target_start = time.time()
            target_cpd_requests = []
            # Write target compounds to target collection
            # Target compounds are written as mine compounds
            print("--------------- Targets ----------------")
            # Insert target compounds
            target_start = time.time()
            # non-parallel insertion
            for comp_dict in self.compounds.values():
                if comp_dict['_id'].startswith('T'):
                    db.insert_mine_compound(comp_dict, target_cpd_requests)     
            print(f"Done with Target Prep--took {time.time() - target_start} seconds.")
            if target_cpd_requests:
                target_start = time.time()
                db.target_compounds.bulk_write(target_cpd_requests, ordered=False)
                print(f"Inserted {len(target_cpd_requests)} Target Compounds in {time.time() - target_start} seconds.")
                del(target_cpd_requests)
                db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                    "Action": "Target Compounds Inserted"})
            else:
                print('No Target Compounds Inserted')
            print(f"----------------------------------------\n")                   

        # Save operators  
        operator_start = time.time()
        if self.operators:
            print("-------------- Operators ---------------")
            # update operator rxn count
            for rxn_dict in self.reactions.values():
                for op in rxn_dict['Operators']:
                    self.operators[op][1]['Reactions_predicted'] += 1
            db.operators.insert_many([op[1] for op in self.operators.values()])        
            db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                    "Action": "Operators Inserted"})  
            print(f"Done with Operators Overall--took {time.time() - operator_start} seconds.")
        print(f"----------------------------------------\n")

        if indexing:
            print("-------------- Indices ---------------")
            index_start = time.time()
            db.build_indexes()
            print(f"Done with Indices--took {time.time() - index_start} seconds.")
            print(f"----------------------------------------\n")

        print("-------------- Overall ---------------")
        print(f"Finished uploading everything in {time.time() - start} sec")
        print(f"----------------------------------------\n")

    def transform(self):

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            n = max(n, 1)           
            for i in range(0, len(lst), n):
                yield lst[i:i + n]

        new_cpds = dict()
        new_rxns = dict()
        compound_smiles = [cpd['SMILES'] for cpd in self.compounds.values()
                            if cpd['Generation'] == self.generation
                            and cpd['Type'] not in ['Coreactant', 'Target Compound']
                            and cpd['Expand'] == True]

        coreactant_dict = {co_key: self.compounds[co_key] for _, co_key in self.coreactants.values()}
        
        self.generation += 1
        for cpd_chunk in chunks(compound_smiles, 3):
            new_cpds_chunk, new_rxns_chunk = _transform_all_helper(cpd_chunk, self.coreactants, coreactant_dict, 
                self.operators, self.generation, self.explicit_h)

            new_cpds.update(new_cpds_chunk)
            for rxn, rxn_dict in new_rxns_chunk.items():
                if rxn in new_rxns:
                    new_rxns[rxn]['Operators'].update(new_rxns[rxn]['Operators'])
                else:
                    new_rxns.update({rxn:rxn_dict})

        # update self.compounds / self.reactions here
        for cpd_id, cpd_dict in new_cpds.items():
            if cpd_id not in self.compounds:
                self.compounds[cpd_id] = cpd_dict
        
        for rxn_id, rxn_dict in new_rxns.items():
            if rxn_id not in self.reactions:
                self.reactions[rxn_id] = rxn_dict
            else:
                self.reactions[rxn_id].update(rxn_dict['Operators'])

            # Update compound tracking
            for prod_id in [cpd_id for _, cpd_id in rxn_dict['Products'] if cpd_id.startswith('C')]:
                if rxn_id not in self.compounds[prod_id]['Product_of']:
                    #TODO make set
                    self.compounds[prod_id]['Product_of'].append(rxn_id)
            
            for reac_id in [cpd_id for _, cpd_id in rxn_dict['Reactants'] if cpd_id.startswith('C')]:
                if rxn_id not in self.compounds[reac_id]['Reactant_in']:
                    self.compounds[reac_id]['Reactant_in'].append(rxn_id)


def _racemization(compound, max_centers=3, carbon_only=True):
    """Enumerates all possible stereoisomers for unassigned chiral centers.

    :param compound: A compound
    :type compound: rdMol object
    :param max_centers: The maximum number of unspecified stereocenters to
        enumerate. Sterioisomers grow 2^n_centers so this cutoff prevents lag
    :type max_centers: int
    :param carbon_only: Only enumerate unspecified carbon centers. (other
        centers are often not tautomeric artifacts)
    :type carbon_only: bool
    :return: list of stereoisomers
    :rtype: list of rdMol objects
    """
    new_comps = []
    # FindMolChiralCenters (rdkit) finds all chiral centers. We get all
    # unassigned centers (represented by '?' in the second element
    # of the function's return parameters).
    unassigned_centers = [c[0] for c in AllChem.FindMolChiralCenters(
        compound, includeUnassigned=True) if c[1] == '?']
    # Get only unassigned centers that are carbon (atomic number of 6) if
    # indicated
    if carbon_only:
        unassigned_centers = list(
            filter(lambda x: compound.GetAtomWithIdx(x).GetAtomicNum() == 6,
                unassigned_centers))
    # Return original compound if no unassigned centers exist (or if above
    # max specified (to prevent lag))
    if not unassigned_centers or len(unassigned_centers) > max_centers:
        return [compound]
    for seq in itertools.product([1, 0], repeat=len(unassigned_centers)):
        for atomid, clockwise in zip(unassigned_centers, seq):
            # Get both cw and ccw chiral centers for each center. Used
            # itertools.product to get all combinations.
            if clockwise:
                compound.GetAtomWithIdx(atomid).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:
                compound.GetAtomWithIdx(atomid).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        # Duplicate C++ object so that we don't get multiple pointers to
        # same object
        new_comps.append(deepcopy(compound))
    return new_comps

def _transform_compound(cpd_smiles, coreactant_mols, coreactant_dict, operators, generation, explicit_h):
    # kekulize

    def _make_half_rxn(mol_list, rules):
        cpds = {}
        cpd_counter = collections.Counter()

        for mol, rule in zip(mol_list, rules):
            if rule == 'Any':
                cpd_dict = _gen_compound(mol)
                # failed compound
                if cpd_dict == None:
                    return None, None
            else:
                cpd_id = coreactant_mols[rule][1]
                cpd_dict = coreactant_dict[cpd_id]            

            cpds[cpd_dict['_id']] = cpd_dict
            cpd_counter.update({cpd_dict['_id']:1})
            
        atom_counts = collections.Counter()
        for cpd_id, cpd_dict in cpds.items():
            for atom_id, atom_count in cpd_dict['atom_count'].items():
                atom_counts[atom_id] += atom_count*cpd_counter[cpd_id]
                
        return [(stoich, cpds[cpd_id]) for cpd_id, stoich in cpd_counter.items()], atom_counts

    def _gen_compound(mol):
        try:
            if explicit_h:
                mol = AllChem.RemoveHs(mol)
            AllChem.SanitizeMol(mol)
        except:
            return None

        mol_smiles = AllChem.MolToSmiles(mol, True)            
        c_id = utils.compound_hash(mol_smiles, 'Predicted')
        if c_id not in local_cpds:
            cpd_dict = {'ID': None, '_id': c_id, 'SMILES': mol_smiles,
                            'Type': 'Predicted',
                            'Generation': generation,
                            'atom_count': utils._getatom_count(mol),
                            'Reactant_in': [], 'Product_of': [],
                            'Expand': True}
        else:
            cpd_dict = local_cpds[c_id]
    
        return cpd_dict

    

    local_cpds = {}
    local_rxns = {}

    mol = AllChem.MolFromSmiles(cpd_smiles)
    mol = AllChem.RemoveHs(mol)

    if not mol:
        if self.errors:
            raise ValueError(f"Unable to parse: {compound_smiles}")
        else:
            print(f"Unable to parse: {compound_smiles}")
            return
    # if self.kekulize:
    #         AllChem.Kekulize(mol, clearAromaticFlags=True)
    # if self.explicit_h:
    #     mol = AllChem.AddHs(mol)

    # Apply reaction rules to prepared compound

    for rule_name, rule in operators.items():
        # Get RDKit Mol objects for reactants
        # TODO: add functionality for bimolecular+ reactions
        reactant_mols = tuple([mol if x == 'Any'
                                else coreactant_mols[x][0]
                                for x in rule[1]['Reactants']])
        # Perform chemical reaction on reactants for each rule
        try:
            product_sets = rule[0].RunReactants(reactant_mols)
        # This error should be addressed in a new version of RDKit
        # TODO: Check this claim
        except RuntimeError:
            print("Runtime ERROR!" + rule_name)
            print(compound_smiles)
            continue
        
        #FIX
        reactants, reactant_atoms = _make_half_rxn(reactant_mols, rule[1]['Reactants'])      
             
        if not reactants:
            continue

        for product_mols in product_sets:
                try:
                    try:
                        #FIX
                        products, product_atoms = _make_half_rxn(product_mols, rule[1]['Products'])
                        if not products:
                            continue
                    except:
                        continue
                    
                    if (reactant_atoms - product_atoms or product_atoms - reactant_atoms):
                        is_atom_balanced = False
                    else:
                        is_atom_balanced = True

                    if is_atom_balanced:
                        for _, cpd_dict in products:
                            if cpd_dict['_id'].startswith('C'):
                                local_cpds.update({cpd_dict['_id']:cpd_dict})
                        #FIX
                        #add reaction here
                        rhash, rxn_text = rxn2hash(reactants, products)
                        if rhash not in local_rxns:
                            local_rxns[rhash] = {'_id': rhash,
                                                 # give stoich and id of reactants/products
                                                 'Reactants': [(s, r['_id']) for s, r in reactants],
                                                 'Products': [(s, p['_id']) for s, p in products],
                                                 'Operators': {rule_name},                                    
                                                 'SMILES_rxn': rxn_text}
                        else:
                            local_rxns[rhash]['Operators'].add(rule_name)                                                                                 
            
                except (ValueError, MemoryError) as e:
                    continue    
    
    return local_cpds,local_rxns

def _transform_all_helper(cpd_list, coreactants, coreactant_dict, operators, generation, explicit_h, **kwargs):
    """
    This function is made to reduce the memory load of parallelization.
    Currently it is believed that the in pickaxe class parallelization will generation
    N copies of the class in memory for each new process, leading to a lot of unecessary memory.

    This function will accept only the necessary information and return the calculated information to be 
    processed later within the class to resolve collisions in cpd_ids and to update product/reactant in.

    This function accepts in a list of cpds (cpd_list) and runs the transformation in parallel of these.
    """
    # process kwargs

    new_cpds_master = {}
    new_rxns_master = {}
    # par loop
    for cpd_smiles in cpd_list:
        new_cpds, new_rxns = _transform_compound(cpd_smiles, coreactants, 
            coreactant_dict, operators, generation, explicit_h)
        # new_cpds as cpd_id:cpd_dict
        # new_rxns as rxn_id:rxn_dict
        new_cpds_master.update(new_cpds)
        # Need to check if reactions already exist to update operators list
        for rxn, rxn_dict in new_rxns.items():
            if rxn in new_rxns_master:
                new_rxns_master[rxn]['Operators'].union(rxn_dict['Operators'])
            else:
                new_rxns_master.update({rxn:rxn_dict})


    return new_cpds_master, new_rxns_master