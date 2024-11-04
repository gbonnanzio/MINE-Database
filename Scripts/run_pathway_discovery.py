import sys
# below path should point to the directory that contains utils directory (and minedatbase if local)
sys.path.append("/projects/p30041/gbf4422/MINE-Database/") 
from utils import utils
import pymongo
import pickle
import json
import os
from minedatabase.pickaxe import Pickaxe
from minedatabase.filters import (SimilarityFilter,SimilaritySamplingFilter)

# file paths in this run file assume that you are running a bash script submission in the scripts directory 

# Weight function for sampling intermediate metabolites by similarity
def weight(score):
    """weight is a function that accepts a similarity score as the sole argument
    and returns a scaled value.
    """
    return score ** 4

if __name__ == '__main__':

### Section 0 - Enter expansion parameters (be sure to check everything here)

    # Set unique expansion ID and describe this expansion
    sub_indx =  0 #int(sys.argv[1])
    exp_ID = 'Reverse_Pyruvate_to_3HPA_' + str(sub_indx) 
    remarks = 'Reverse_Pyruvate_to_3HPA_'+ str(sub_indx)

    # Enter precursor and target names as well as SMILES
    precursor_name = '3HPA' + str(sub_indx)
    precursor_list = ['OCC(C)(C)C(O)=O']
    precursor_smiles = precursor_list[sub_indx]

    # if no target put None
    target_name = 'Pyruvate'
    target_smiles = 'OC(C(C)=O)=O'

    # Set number of generations and processes for biochemical network expansion
    num_generations = 2
    num_processes = 40

    # Set similarity sampling filter and number of compounds to sample per generation
    similarity_sample = False
    num_samples = 1000

    # Set similarity thresholds and values to cutoff compounds at per generation
    similarity_filter = False
    increasing_similarity = False # leave on False, then compounds don't need to strictly increase in similarity
    similarity_threshold = [0.0, 0.0, 0.0]

    # Pick cofactors to expand with
    #coreactant_list = '../data/coreactants_and_rules/all_cofactors.tsv'
    coreactant_list = '../minedatabase/data/A_plus_B_rules/COFACTORS_GB_4784_rules_543_cofactors_DMB.tsv'

    # Decide if saving results locally or on mongo
    write_mongo = False
    write_local = True

    mongo_conn = 'mongodb://ymc1840:yash_is_cool@minedatabase.ci.northwestern.edu:27017'
    local_dir = '../minedatabase/data/pickaxe_runs/'

    # Decide if you want to visualize the reaction network as an interactive bokeh plot
    # (will take long time if >2 generations)
    generate_interactive_plot = False

### Section 1 - Initialize Pickaxe

    # canonicalize precursor and target SMILES then write to tsv
    precursor_smiles = utils.canonicalize_smiles(precursor_smiles)
    precursor_filepath = utils.write_cpds_to_tsv(precursor_name,precursor_smiles)

    if target_smiles is not None and target_name is not None:
        target_smiles = utils.canonicalize_smiles(target_smiles)
        target_filepath = utils.write_cpds_to_tsv(target_name,target_smiles)

    # pick relevant rules for expansion ('generalized', 'intermediate', or 'intermediate_non_dimerization')
    # or specify a custom list for A+B expansions by specifying that rule's path
    rules_type = '../minedatabase/data/A_plus_B_rules/RULES_GB_4784_rules_543_cofactors_DMB.tsv'

    rules_range = None # e.g. 100 selects the top 100 rules, None means use all rules
    specific_rule = None # e.g. 'rule0004' or 'rule0004_03', None means use all rules

    # if we specified a custom rule set (like for A+B expansions)
    if '.tsv' in rules_type: 
        rule_list = rules_type
    else:
        rule_list = utils.pick_rules(rules_type = rules_type,
                                 rules_range = rules_range,
                                 specific_rule = specific_rule)

    # initialize pickaxe object
    pk = Pickaxe(coreactant_list = coreactant_list, rule_list = rule_list)
    pk.load_compound_set(compound_file=precursor_filepath) # load input compound in Pickaxe
    if target_smiles is not None and target_name is not None:
        pk.load_targets(target_filepath) # load target compound in Pickaxe

    ## incorporate filters into pickaxe

    # Similarity filter
    sample_fingerprint_method = "Morgan"
    cutoff_fingerprint_method = "Morgan"
    cutoff_fingerprint_args = {"radius": 2}
    cutoff_similarity_method = "Tanimoto"

    crit_similarity = taniFilter = SimilarityFilter(
                crit_similarity = similarity_threshold,
                increasing_similarity=increasing_similarity,
                fingerprint_method=sample_fingerprint_method,
                fingerprint_args=cutoff_fingerprint_args,
                similarity_method=cutoff_similarity_method)

    pk.filters.append(crit_similarity)

    # Similarity sampling filter
    sample_size = num_samples
    sample_fingerprint_method = "Morgan"
    sample_fingerprint_args = {"radius": 2}
    sample_similarity_method = "Tanimoto"

    if similarity_sample:
        taniSampleFilter = SimilaritySamplingFilter(
            sample_size=sample_size,
            weight=weight,
            fingerprint_method=sample_fingerprint_method,
            fingerprint_args=sample_fingerprint_args,
            similarity_method=sample_similarity_method)
        pk.filters.append(taniSampleFilter)

    # Save all details about this expansion locally or on MongoDB
    expansion_details = {"expansion_ID": exp_ID,
                         "precursor_name": precursor_name,
                         "precursor_SMILES": precursor_smiles,
                         "target_name": target_name,
                         "target_SMILES": target_smiles,
                         "num_generations": num_generations,
                         "num_processes": num_processes,
                         "rules_type": rules_type,
                         "rules_range": rules_range,
                         "specific_rule": specific_rule,
                         "similarity_sample": similarity_sample,
                         "similarity_sampling_size": num_samples,
                         "similarity_file": similarity_filter,
                         "similarity_thresholds": similarity_threshold,
                         "remarks": remarks}

    # If user chooses to write expansion details locally
    if write_local:

        try:
            # Create a new folder with the same name as the expansion ID
            os.mkdir(f"../minedatabase/data/pickaxe_runs/{exp_ID}")

        except FileExistsError:
            pass

            # Save expansion details locally
        with open(f"../minedatabase/data/pickaxe_runs/{exp_ID}/{exp_ID}_expansion_details.json", "w") as outfile:
            json.dump(expansion_details, outfile)

    # If user chooses to write expansion details on mongo
    if write_mongo:

        # Connect to MongoDB, then create a db and collection
        mongo_client = pymongo.MongoClient(mongo_conn)
        this_expansion_db = mongo_client[exp_ID]
        exp_details_col = this_expansion_db['expansion_details']
        docs_alr_present = exp_details_col.find({})

        # Check if there are any documents already present in this collection
        i = 0
        for doc in docs_alr_present:
            i += 1

        if i >= 1:
            raise Exception("Please define a different expansion ID. This one already exists")

        else:
            exp_details_col.insert_one(expansion_details)

### Section 2 - Perform Pickaxe expansion

    pk.transform_all(generations=num_generations,processes=num_processes)
    pk.assign_ids()

    if write_local:
        # then store all compounds and reactions associated with the pickaxe object
        pk.pickle_pickaxe(f'../minedatabase/data/pickaxe_runs/{exp_ID}/{exp_ID}.pkl')
        #with open(f'../data/pickaxe_runs/{exp_ID}/{exp_ID}.pkl', 'wb') as file:  # 'wb' for write binary
        #    pickle.dump(pk, file)

    if write_mongo:

        mongo_client = pymongo.MongoClient(mongo_conn)
        this_expansion_db = mongo_client[exp_ID]
        compounds_col = this_expansion_db['compounds']
        reactions_col = this_expansion_db['reactions']

        for compound in pk.compounds:
            compounds_col.insert_one(pk.compounds[compound])

        for reaction in pk.reactions:
            rxn_dict = pk.reactions[reaction]
            rxn_entry = {'_id':rxn_dict['_id'],
                         'Reactants': rxn_dict['Reactants'],
                         'Products': rxn_dict['Products'],
                         'Operators': list(rxn_dict['Operators']),
                         'Reaction (SMILES)' : rxn_dict['SMILES_rxn'],
                         'Reaction (IDs)': rxn_dict['ID_rxn'],
                         'Reaction ID': rxn_dict['ID']}

            reactions_col.insert_one(rxn_entry)
                                      
print('\n----- Finished --------')