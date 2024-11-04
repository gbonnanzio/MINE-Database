from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import matplotlib.pyplot as plt
from IPython.display import Image, display
import pandas as pd
import ast

def draw_molecule(smiles,size_tuple):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        # Generate an image of the molecule
        img = Draw.MolToImage(mol, size_tuple)
        # Display the image
        display(img)
    else:
        print("Invalid SMILES string.")

_REACTIONS = None
def neutralise_charges(
    mol: Chem.rdchem.Mol, reactions=None
) -> Chem.rdchem.Mol:
    """Neutralize all charges in an rdkit mol.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule to neutralize.
    reactions : list, optional
        patterns to neutralize, by default None.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        Neutralized molecule.
    """

    def _initialise_neutralisation_reactions():
        patts = (
            # Imidazoles
            ("[n+;H]", "n"),
            # Amines
            ("[N+;!H0]", "N"),
            # Carboxylic acids and alcohols
            ("[$([O-]);!$([O-][#7])]", "O"),
            # Thiols
            ("[S-;X1]", "S"),
            # Sulfonamides
            ("[$([N-;X2]S(=O)=O)]", "N"),
            # Enamines
            ("[$([N-;X2][C,N]=C)]", "N"),
            # Tetrazoles
            ("[n-]", "[nH]"),
            # Sulfoxides
            ("[$([S-]=O)]", "S"),
            # Amides
            ("[$([N-]C=O)]", "N"),
        )
        return [
            (AllChem.MolFromSmarts(x), AllChem.MolFromSmiles(y, False))
            for x, y in patts
        ]

    global _REACTIONS  # pylint: disable=global-statement
    if reactions is None:
        if _REACTIONS is None:
            _REACTIONS = _initialise_neutralisation_reactions()
        reactions = _REACTIONS
    for (reactant, product) in reactions:
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol

def clean_smiles(smiles,remove_stereo,neutralize):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if remove_stereo == True:
            Chem.rdmolops.RemoveStereochemistry(mol)
            # Remove stereochemistry information        
        if neutralize == True:
            # Neutralize charges on atoms
            mol = neutralise_charges(mol)
            cpd_smiles = Chem.MolToSmiles(mol)
            mol = Chem.MolFromSmiles(cpd_smiles)        
        # Generate a new SMILES string without stereochemistry and with neutral atoms
        cpd_smiles = Chem.MolToSmiles(mol)
        if mol == None:
            return None
        else:
            return cpd_smiles
    except:
        return None

def draw_reaction(reactant_smiles,product_smiles,size):
    # Convert the smiles to RDKit molecules
    reactants = [Chem.MolFromSmiles(smiles) for smiles in reactant_smiles]
    products = [Chem.MolFromSmiles(smiles) for smiles in product_smiles]

    # Generate the reaction SMARTS
    reactant_smarts = '.'.join([Chem.MolToSmarts(reactant) for reactant in reactants])
    product_smarts = '.'.join([Chem.MolToSmarts(product) for product in products])
    reaction_smarts = f'{reactant_smarts}>>{product_smarts}'

    # Generate the reaction
    reaction = Chem.rdChemReactions.ReactionFromSmarts(reaction_smarts)

    # Generate the reaction image
    reaction_image = Draw.ReactionToImage(reaction)
    #display(reaction_image)
    plt.figure(figsize=size)
    plt.imshow(reaction_image)
    plt.axis('off')
    plt.show()

def get_SMILES_from_alias(SEED_df,alias):
    """
    Get SMILES string and ID of a compound from its name
    """

    count = 0
    idx = []
    for i in range(0,SEED_df.shape[0]):
        alias_strings = SEED_df['aliases'][i]
        aliases_list = ast.literal_eval(alias_strings) #alias_strings.strip("[' ']").split(',')
        for curr_alias in aliases_list:
            if alias.lower() == curr_alias.lower().strip():
                idx.append(i)
                count += 1
                break
        if len(idx) > 0:
            break
    
    if count != 0: # if compound is located
        cpd_smiles = SEED_df['SMILES'].iloc[idx[0]]
    
        cpd_id =SEED_df["ID"].iloc[idx[0]]
    
    else: # if compound is not located:
        cpd_smiles = None
        cpd_id = None
        
    return cpd_smiles,cpd_id,count