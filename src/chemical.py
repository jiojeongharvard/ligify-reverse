import requests
from rdkit import Chem
from rdkit.Chem import Draw
import re
from libchebipy import ChebiEntity


def get_smiles_from_chebi(chebi_id):
    smiles = ChebiEntity("CHEBI:15372").get_smiles()
    return smiles


def smiles_to_image(smiles, img_file='chemical_structure.png'):
    """Generate an image from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        img.save(img_file)
        return img
    else:
        print("Failed to generate image.")
        return None