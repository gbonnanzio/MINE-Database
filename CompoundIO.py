#!/usr/bin/env python
"""Databases.py: Functions to load MINE databases from and dump compounds into common cheminformatics formats"""

__author__ = 'JGJeffryes'

from pymongo import DESCENDING
from rdkit.Chem import AllChem
from Databases import MINE
import Utils
import sys
import hashlib


def export_sdf(mine_db, target):
    """
    Exports compounds from the database as an MDL SDF file
    :param mine_db: a MINE object
    :param path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to 100000)
    :return:
    """
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as an SDF file") %(mine_db.compounds.count(), mine_db.name)
    w = AllChem.SDWriter(target)
    for compound in mine_db.compounds.find():
        mol = AllChem.MolFromSmiles(compound['SMILES'], True, {'CoA': '*', 'R': "*"})
        mol.SetProp('_id', compound['_id'])
        mol.SetProp('Generation', compound['Generation'])
        if 'Reactant_in' in compound:
            mol.SetProp('Reactant_in', compound['Reactant_in'])
        if 'Reactant_in' in compound:
            mol.SetProp('Product_of', compound['Product_of'])
        w.write(mol)
    w.close()

def import_sdf(mine_db, target, id_db='UniversalMINE'):
    """
    Imports a SDF file as a MINE database
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the SDF file to be loaded
    :param id_db: string, name of the db which coordinates MINE_ids
    :return:
    """
    last_cid = None
    if id_db:
        uni_mine = MINE(id_db)
        last_cid = [int(x['MINE_id']) for x in uni_mine.compounds.find().sort([("MINE_id", DESCENDING)]).limit(1)][0]

    sdf_gen = AllChem.SDMolSupplier(target)
    for mol in sdf_gen:
        compound = {}
        for key in mol.GetPropNames():
            compound[key] = mol.GetProp(key)
        compound['SMILES'] = AllChem.MolToSmiles(mol)
        compound['Inchikey'] = AllChem.MolToInchikey(mol)
        compound['Inchi'] = AllChem.MolToInchi(mol)
        compound['Mass'] = AllChem.CalcExactMolWt(mol)
        compound['Formula'] = AllChem.CalcMolFormula(mol)
        compound['Charge'] = AllChem.GetFormalCharge(mol)
        compound['MACCS'] = [i for i, bit in enumerate(AllChem.GetMACCSKeysFingerprint(mol)) if bit]
        compound['len_MACCS'] = len(compound['MACCS'])
        compound['RDKit'] = [i for i, bit in enumerate(AllChem.RDKFingerprint(mol)) if bit]
        compound['len_RDKit'] = len(compound['RDKit'])

        comphash = hashlib.sha1(compound['SMILES']).hexdigest()
        if '_id' in compound:
            if "X" in compound['_id']:
                compound = Utils.fix_rxn_pointers(mine_db, 'X' + comphash, compound)
            else:
                compound = Utils.fix_rxn_pointers(mine_db, 'C' + comphash, compound)
        else:
            compound['_id'] = 'C' + comphash

        compound['DB_links'] = {}
        #TODO: insert DB_links

        if id_db:
            mine_comp = uni_mine.compounds.find_one({"_id": compound['_id']})
            if mine_comp:
                compound['MINE_id'] = mine_comp['MINE_id']
            elif compound['_id'][0] == 'C':
                last_cid += 1
                compound['MINE_id'] = last_cid
                uni_mine.compounds.insert(Utils.convert_sets_to_lists(compound))

        mine_db.compounds.save(Utils.convert_sets_to_lists(compound))


if __name__ == '__main__':
    task = sys.argv[1]
    db_name = sys.argv[2]
    path = sys.argv[3]
    database = MINE(db_name)
    if task == 'export-sdf':
        export_sdf(database, path)
    if task == 'import-sdf':
        import_sdf(database, path)
