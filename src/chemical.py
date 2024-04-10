import requests
from rdkit import Chem
from rdkit.Chem import Draw
import re

def get_smiles_from_chebi(chebi_id):
    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
    response = requests.get(url)
    if response.status_code == 200:
        html_content = response.text
        # Adjusted regular expression to match the provided HTML snippet
        match = re.search(r'<td class="chebiDataHeader".*?>SMILES</td>\s*<td>(.*?)</td>', html_content, re.DOTALL)
        if match:
            result = match.group(1).strip()
            result = result.replace('<img alt="" src="images/invisible_space.gif" width="0" height="0" border="0"/>', '')
            return result
    else:
        print("Error fetching the data")
        return None


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