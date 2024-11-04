import csv
import pandas as pd
from rdkit import Chem, RDLogger
import networkx as nx
import re
import pymongo
import itertools
import json
import numpy as np

# Silence non-critical RDKit warnings to minimize unnecessary outputs
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def canonicalize_smiles(smi: str):
    """
    Canonicalize SMILES string with RDKit
    :param: Input SMILES string, assume not canonical
    :return: Canonicalized SMILES string
    """

    uncanon_smi = smi  # assume input SMILES string is not canonical

    try:  # try canonicalizing first
        canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(uncanon_smi))

    except:  # leave SMILE as is if it cannot be canonicalized
        canon_smi = uncanon_smi

    return canon_smi


def search_SMILES_in_biological_databases(SMILES: str, biological_SMILES: set):
    """
    Search SMILES of a compound against combined BRENDA, KEGG, and Metacyc databases
    @param SMILES: SMILES of queried compound
    @param biological_SMILES: list of SMILES present in BRENDA, KEGG, and Metacyc databases
    @return: True if compound is present in biological databases
    """
    canonicalized_SMILES = canonicalize_smiles(SMILES)  # canonicalize queried SMILES
    return (
        canonicalized_SMILES in biological_SMILES
    )  # True if compound is present in biological databases


def write_cpds_to_tsv(cpd_name: str, cpd_smi: str):
    """
    Write a compound name and SMILES to a .tsv extension file for use in Pickaxe
    :param cpd_name: compound name
    :param cpd_smi: compound SMILES (does not need to be canonicalized)
    :return: filepath of this compound
    """

    cpd_filepath = "../minedatabase/data/pickaxe_runs/" + cpd_name + ".tsv"

    with open(cpd_filepath, "wt", encoding="utf-8") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(["id", "smiles"])
        tsv_writer.writerow([cpd_name, cpd_smi])

    return cpd_filepath


def pick_rules(rules_type: str, rules_range: int = None, specific_rule: str = None):
    path = '../minedatabase/data/original_rules/'
    # If user wants to expand with all 1224 generalized operators
    if rules_type == "generalized" and rules_range is None and specific_rule is None:
        rules_path = path + 'JN1224MIN_rules.tsv'
        return rules_path

    # If user wants to expand with all 3604 intermediate operators
    if rules_type == "intermediate" and rules_range is None and specific_rule is None:
        rules_path = path + "JN3604IMT_rules.tsv"
        return rules_path

    # If user wants to expand with intermediate operators but exclude dimerization operators
    if rules_type == 'intermediate_non_dimerization' and rules_range is None and specific_rule is None:
        rules_path = path + "non_dimerization_imt_rules.tsv"
        return rules_path

    # If user wants to expand with the top X generalized operators starting from rule0001
    if (
        rules_type == "generalized"
        and rules_range is not None
        and specific_rule is None
    ):
        # Read in generalized rules first as dataframe
        rules_path = path + "JN1224MIN_rules.tsv"
        rules_df = pd.read_csv(rules_path, delimiter="\t")

        # Then slice the relevant top X rules starting from rule0001
        input_rules_df = rules_df.iloc[0:rules_range, :].set_index("Name")

        # Save this new dataframe and return its filepath
        new_rules_path = path + "input_rules.tsv"
        input_rules_df.to_csv(new_rules_path, sep="\t")

        return new_rules_path

    # If user wants to expand with a specific generalized operator
    if (
        rules_type == "generalized"
        and rules_range is None
        and specific_rule is not None
    ):
        # Read in generalized rules first as dataframe
        rules_path = path + "JN1224MIN_rules.tsv"
        rules_df = pd.read_csv(rules_path, delimiter="\t")

        # Then index out the relevant rule
        input_rules_df = rules_df[rules_df["Name"] == specific_rule].set_index("Name")

        # Save this new dataframe and return its filepath
        new_rules_path = path + "input_rules.tsv"
        input_rules_df.to_csv(new_rules_path, sep="\t")

        return new_rules_path

    # If user wants to expand with top X intermediate operators starting from rule0001_001
    if (
        rules_type == "intermediate"
        and rules_range is not None
        and specific_rule is None
    ):
        # Read in generalized rules first as dataframe
        rules_path = path + "JN3604IMT_rules.tsv"
        rules_df = pd.read_csv(rules_path, delimiter="\t")

        # Then slice the relevant top X rules starting from rule0001
        input_rules_df = rules_df.iloc[0:rules_range, :].set_index("Name")

        # Save this new dataframe and return its filepath
        new_rules_path = path + "input_rules.tsv"
        input_rules_df.to_csv(new_rules_path, sep="\t")

        return new_rules_path

        # If user wants to expand with a specific generalized operator

    # If user wants to expand with a specific intermediate operator
    if (
        rules_type == "intermediate"
        and rules_range is None
        and specific_rule is not None
    ):
        # Read in generalized rules first as dataframe
        rules_path = path + "JN3604IMT_rules.tsv"
        rules_df = pd.read_csv(rules_path, delimiter="\t")

        # Then index out the relevant rule
        input_rules_df = rules_df[rules_df["Name"] == specific_rule].set_index("Name")

        # Save this new dataframe and return its filepath
        new_rules_path = path + "input_rules.tsv"
        input_rules_df.to_csv(new_rules_path, sep="\t")

        return new_rules_path

    return None


def create_compounds_df(pk: any):
    """
    Compile compounds generated by Pickaxe (includes precursor and cofactors)
    :param pk: Pickaxe object
    :return: Pandas Dataframe of compounds
    """
    pk_cpd_keys = [key for key in pk.compounds.keys()]

    cpd_id_list = [pk.compounds[key]["ID"] for key in pk_cpd_keys]
    cpd_type_list = [pk.compounds[key]["Type"] for key in pk_cpd_keys]
    cpd_gen_list = [pk.compounds[key]["Generation"] for key in pk_cpd_keys]
    cpd_formula_list = [pk.compounds[key]["Formula"] for key in pk_cpd_keys]
    cpd_smiles_list = [pk.compounds[key]["SMILES"] for key in pk_cpd_keys]
    cpd_reactant_list = [
        pk.compounds[key]["Reactant_in"]
        if "Reactant_in" in pk.compounds[key].keys()
        else ""
        for key in pk_cpd_keys
    ]
    cpd_product_list = [
        pk.compounds[key]["Product_of"]
        if "Product_of" in pk.compounds[key].keys()
        else ""
        for key in pk_cpd_keys
    ]

    # Combine all lists into a Pandas Dataframe
    compounds_df = pd.DataFrame(
        {
            "ID": cpd_id_list,
            "Type": cpd_type_list,
            "Generation": cpd_gen_list,
            "Formula": cpd_formula_list,
            "SMILES": cpd_smiles_list,
            "Reactant_in": cpd_reactant_list,
            "Product_in": cpd_product_list,
        }
    )

    return compounds_df


def save_pk_rxns_locally(pk: any, exp_ID: str):
    all_rxns = {}

    for reaction in pk.reactions:
        rxn_dict = pk.reactions[reaction]

        new_dict = {
            rxn_dict["_id"]: {
                "_id": rxn_dict["_id"],
                "Reactants": rxn_dict["Reactants"],
                "Products": rxn_dict["Products"],
                "Operators": str(rxn_dict["Operators"]),  # convert set to str
                "SMILES_rxn": rxn_dict["SMILES_rxn"],
                "ID_rxn": rxn_dict["ID_rxn"],
                "ID": rxn_dict["ID"],
            }
        }

        all_rxns.update(new_dict)

    with open(f"../minedatabase/data/pickaxe_runs/{exp_ID}/{exp_ID}_reactions.json", "w") as outfile:
        json.dump(all_rxns, outfile)


def save_pk_cpds_locally(pk: any, exp_ID: str):
    with open(f"../minedatabase/data/pickaxe_runs/{exp_ID}/{exp_ID}_compounds.json", "w") as outfile:
        json.dump(pk.compounds, outfile)


def create_graph(all_rxn_strs_in_cpd_ids: list, precursor_name: str):
    # initialize graph
    G = nx.DiGraph()

    for i, rxn in enumerate(all_rxn_strs_in_cpd_ids):
        rxn_lhs = rxn.split(" => ")[0]
        rxn_rhs = rxn.split(" => ")[1]

        pattern = r"^(\(\d\)\s)?(\S+)(\[c0\])$"

        substrates_list = [
            re.sub(pattern, r"\2", substrate) for substrate in rxn_lhs.split(" + ")
        ]
        substrates_list = [
            substrate
            for substrate in substrates_list
            if "pkc" in substrate or precursor_name == substrate
        ]

        products_list = [
            re.sub(pattern, r"\2", product) for product in rxn_rhs.split(" + ")
        ]
        products_list = [product for product in products_list if "pkc" in product]

        for substrate in substrates_list:
            for product in products_list:
                G.add_edge(substrate, product)

    return G


def get_sequences_from_graph(
    G, compounds_df, precursor_smiles: str, target_smiles: str, num_generations: int
):
    try:
        source_node = list(
            compounds_df[compounds_df["SMILES"] == precursor_smiles]["ID"]
        )[0]
        target_node = list(compounds_df[compounds_df["SMILES"] == target_smiles]["ID"])[
            0
        ]

        paths = nx.all_simple_paths(G, source_node, target_node, cutoff=num_generations)

        pathway_list = list(paths)
        return pathway_list

    except IndexError:
        print(
            "No pathways were found. Check the compounds file for other potential targets"
        )
        exit()  # terminate pathway discovery altogether if no pathways exist


def get_seq_info(seq: list, biological_compounds: set, compounds_df):
    # count number of intermediates in pathway (excludes precursor and target)
    num_intermediates_total = len(seq[1:-1])

    # initialize a counter to track if intermediates were observed in BRENDA, KEGG, or METACYC
    num_known_intermediates = 0

    # initialize an empty list to track which intermediates were observed in BRENDA, KEGG, or METACYC
    known_intermediate_or_not = []

    for intermediate in seq[1:-1]:
        cpd_SMILES = list(compounds_df[compounds_df["ID"] == intermediate]["SMILES"])[0]
        cpd_reported_in_dbs = search_SMILES_in_biological_databases(
            cpd_SMILES, biological_compounds
        )

        if cpd_reported_in_dbs == True:
            num_known_intermediates += 1
            known_intermediate_or_not.append(cpd_reported_in_dbs)

        else:
            known_intermediate_or_not.append(cpd_reported_in_dbs)

    try:
        proportion_known_intermediates = (
            num_known_intermediates / num_intermediates_total
        )

    except ZeroDivisionError:
        proportion_known_intermediates = None

    return (
        num_intermediates_total,
        num_known_intermediates,
        proportion_known_intermediates,
        known_intermediate_or_not,
    )


def get_pathways_from_graph(
    write_mongo: bool,
    write_local: bool,
    mongo_conn,
    exp_ID: str,
    sequences: list,
    biological_compounds: set,
    compounds_df: any,
    pk: any):

    seq_num = 0
    pathway_num = 0
    all_pathways = {}

    if write_mongo:
        mongo_client = pymongo.MongoClient(mongo_conn)
        this_expansion_db = mongo_client[exp_ID]
        pathways_col = this_expansion_db["pathways"]

    ### Sequence level
    for seq in sequences:
        seq_num += 1
        all_rxns_in_this_seq = []

        # Extract compounds in a pairwise fashion for this sequence
        for i in range(len(seq) - 1):
            pair = seq[i : i + 2]

            substrate_ID = pair[0]  # eg: pkc0000035
            product_ID = pair[1]  # eg: pkc0000059

            # Get reactions that this substrate participates it (index 0 to get values from series object)
            substrate_is_reactant_in = list(
                compounds_df[compounds_df["ID"] == substrate_ID]["Reactant_in"]
            )[0]

            # Get product that this product is formed in (index 0 to get values from series object)
            product_is_product_in = list(
                compounds_df[compounds_df["ID"] == product_ID]["Product_in"]
            )[0]

            # Quick test to ensure all reaction hash keys are unique
            assert len(substrate_is_reactant_in) == len(set(substrate_is_reactant_in))
            assert len(product_is_product_in) == len(set(product_is_product_in))

            # Get all reactions between this substrate and product
            common_rxns = list(
                set(substrate_is_reactant_in).intersection(product_is_product_in)
            )

            # store these reactions
            all_rxns_in_this_seq.append(common_rxns)

        ### Pathway level
        # Get all combinations of reactions between these metabolites
        all_pathways_in_this_seq = list(itertools.product(*all_rxns_in_this_seq))

        for pathway in all_pathways_in_this_seq:
            # update pathway number
            pathway_num += 1

            # initialize empty lists to track info for this pathway
            pathway_rxn_hashes = []
            pathway_rxns_in_SMILES = []
            pathway_rxns_in_cpd_IDs = []
            pathway_rxn_rules = []

            for rxn_hash in pathway:
                # extract pickaxe reaction dict for this reaction in the pathway
                pk_rxn_dict = pk.reactions[rxn_hash]

                # extract the reaction rule for this reaction in the pathway
                rxn_rules = list(pk_rxn_dict["Operators"])

                if len(rxn_rules) == 1:
                    rxn_rules = rxn_rules[0]

                else:
                    rxn_rules = rxn_rules[0] + ";" + rxn_rules[1]

                # extract the reaction string in terms of compound SMILES for this reaction in the pathway
                rxn_str_in_SMILES = pk_rxn_dict["SMILES_rxn"]

                # extract the reaction string in terms of compound IDs for this reaction in the pathway
                rxn_str_in_cpd_IDs = pk_rxn_dict["ID_rxn"]

                # update lists for this pathway
                pathway_rxn_hashes.append(rxn_hash)
                pathway_rxns_in_SMILES.append(rxn_str_in_SMILES)
                pathway_rxns_in_cpd_IDs.append(rxn_str_in_cpd_IDs)
                pathway_rxn_rules.append(rxn_rules)

            # Check how many intermediates in sequence were reported in BRENDA, KEGG, or METACYC
            (
                num_intermediates_total,
                num_known_intermediates,
                proportion_known_intermediates,
                known_intermediates_or_not,
            ) = get_seq_info(seq, biological_compounds, compounds_df)

            # Extract compounds in sequence
            all_seq_cpd_IDs = []
            all_seq_cpd_SMILES = []

            for cpd_ID in seq:
                cpd_SMILES = list(compounds_df[compounds_df["ID"] == cpd_ID]["SMILES"])[
                    0
                ]
                all_seq_cpd_IDs.append(cpd_ID)
                all_seq_cpd_SMILES.append(cpd_SMILES)

            # # Extract compounds in a pairwise fashion for this sequence again to calculate feasibility scores
            #
            # all_feasibility_scores = []

            for i in range(len(seq) - 1):
                pair = seq[i : i + 2]

                substrate_ID = pair[0]
                product_ID = pair[1]

                substrate_SMILES = list(
                    compounds_df[compounds_df["ID"] == substrate_ID]["SMILES"]
                )[0]
                product_SMILES = list(
                    compounds_df[compounds_df["ID"] == product_ID]["SMILES"]
                )[0]

                rxn_rule = list(pathway_rxn_rules)[i]

                if "_" in rxn_rule:
                    gen_rxn_rule = rxn_rule[0:8]

                else:
                    gen_rxn_rule = rxn_rule

                # feasibility_score = PX.predict_feasibility(
                #     substrate_SMILES, product_SMILES, gen_rxn_rule, return_proba=True
                # )

                # all_feasibility_scores.append(feasibility_score)

            #     # convert feasibility scores to strings instead of floats
            # all_feasibility_scores_strs = [
            #     str(feasibility) for feasibility in all_feasibility_scores
            # ]

            pathway_dict = {
                "pathway_num": pathway_num,
                "sequence_num": seq_num,
                "sequence (compound IDs)": all_seq_cpd_IDs,
                "sequence (SMILES)": all_seq_cpd_SMILES,
                "reactions (SMILES)": pathway_rxns_in_SMILES,
                "reactions (compound IDs)": pathway_rxns_in_cpd_IDs,
                "reaction rules": pathway_rxn_rules,
                "reactions (hash keys)": pathway_rxn_hashes,
                "num_intermediates": num_intermediates_total,
                "num_known_intermediates": num_known_intermediates,
                "proportion_known_intermediates": proportion_known_intermediates,
                "type_known_intermediates": known_intermediates_or_not}


            if write_mongo:
                pathways_col.insert_one(pathway_dict)

            if write_local:
                all_pathways.update({f"Pathway {pathway_num}": pathway_dict})

    if write_local:
        with open(
            f"../minedatabase/data/pickaxe_runs/{exp_ID}/{exp_ID}_pathways.json", "w"
        ) as outfile:
            json.dump(all_pathways, outfile)

def get_pathways_from_loaded_graph(
    sequences: list,
    biological_compounds: set,
    compounds_df: any,
    pk: any):

    seq_num = 0
    pathway_num = 0
    all_pathways = {}

    ### Sequence level
    for seq in sequences:
        seq_num += 1
        all_rxns_in_this_seq = []

        # Extract compounds in a pairwise fashion for this sequence
        for i in range(len(seq) - 1):
            pair = seq[i : i + 2]

            substrate_ID = pair[0]  # eg: pkc0000035
            product_ID = pair[1]  # eg: pkc0000059

            # Get reactions that this substrate participates it (index 0 to get values from series object)
            substrate_is_reactant_in = list(compounds_df[compounds_df["ID"] == substrate_ID]["Reactant_in"])[0]

            # Get product that this product is formed in (index 0 to get values from series object)
            product_is_product_in = list(compounds_df[compounds_df["ID"] == product_ID]["Product_in"])[0]

            # Quick test to ensure all reaction hash keys are unique
            assert len(substrate_is_reactant_in) == len(set(substrate_is_reactant_in))
            assert len(product_is_product_in) == len(set(product_is_product_in))

            # Get all reactions between this substrate and product
            common_rxns = list(set(substrate_is_reactant_in).intersection(product_is_product_in))

            # store these reactions
            all_rxns_in_this_seq.append(common_rxns)

        ### Pathway level
        # Get all combinations of reactions between these metabolites
        all_pathways_in_this_seq = list(itertools.product(*all_rxns_in_this_seq))

        for pathway in all_pathways_in_this_seq:
            # update pathway number
            pathway_num += 1

            # initialize empty lists to track info for this pathway
            pathway_rxn_hashes = []
            pathway_rxns_in_SMILES = []
            pathway_rxns_in_cpd_IDs = []
            pathway_rxn_rules = []

            for rxn_hash in pathway:
                # extract pickaxe reaction dict for this reaction in the pathway
                pk_rxn_dict = pk.reactions[rxn_hash]

                # extract the reaction rule for this reaction in the pathway
                rxn_rules = list(pk_rxn_dict["Operators"])

                if len(rxn_rules) == 1:
                    rxn_rules = rxn_rules[0]

                else:
                    rxn_rules = rxn_rules[0] + ";" + rxn_rules[1]

                # extract the reaction string in terms of compound SMILES for this reaction in the pathway
                rxn_str_in_SMILES = pk_rxn_dict["SMILES_rxn"]

                # extract the reaction string in terms of compound IDs for this reaction in the pathway
                rxn_str_in_cpd_IDs = pk_rxn_dict["ID_rxn"]

                # update lists for this pathway
                pathway_rxn_hashes.append(rxn_hash)
                pathway_rxns_in_SMILES.append(rxn_str_in_SMILES)
                pathway_rxns_in_cpd_IDs.append(rxn_str_in_cpd_IDs)
                pathway_rxn_rules.append(rxn_rules)

            # Check how many intermediates in sequence were reported in BRENDA, KEGG, or METACYC
            (
                num_intermediates_total,
                num_known_intermediates,
                proportion_known_intermediates,
                known_intermediates_or_not,
            ) = get_seq_info(seq, biological_compounds, compounds_df)

            # Extract compounds in sequence
            all_seq_cpd_IDs = []
            all_seq_cpd_SMILES = []

            for cpd_ID in seq:
                cpd_SMILES = list(compounds_df[compounds_df["ID"] == cpd_ID]["SMILES"])[
                    0
                ]
                all_seq_cpd_IDs.append(cpd_ID)
                all_seq_cpd_SMILES.append(cpd_SMILES)

            for i in range(len(seq) - 1):
                pair = seq[i : i + 2]

                substrate_ID = pair[0]
                product_ID = pair[1]

                substrate_SMILES = list(
                    compounds_df[compounds_df["ID"] == substrate_ID]["SMILES"]
                )[0]
                product_SMILES = list(
                    compounds_df[compounds_df["ID"] == product_ID]["SMILES"]
                )[0]

                rxn_rule = list(pathway_rxn_rules)[i]

                if "_" in rxn_rule:
                    gen_rxn_rule = rxn_rule[0:8]

                else:
                    gen_rxn_rule = rxn_rule

            pathway_dict = {
                "pathway_num": pathway_num,
                "sequence_num": seq_num,
                "sequence (compound IDs)": all_seq_cpd_IDs,
                "sequence (SMILES)": all_seq_cpd_SMILES,
                "reactions (SMILES)": pathway_rxns_in_SMILES,
                "reactions (compound IDs)": pathway_rxns_in_cpd_IDs,
                "reaction rules": pathway_rxn_rules,
                "reactions (hash keys)": pathway_rxn_hashes,
                "num_intermediates": num_intermediates_total,
                "num_known_intermediates": num_known_intermediates,
                "proportion_known_intermediates": proportion_known_intermediates,
                "type_known_intermediates": known_intermediates_or_not}


            all_pathways.update({f"Pathway {pathway_num}": pathway_dict})

    return all_pathways