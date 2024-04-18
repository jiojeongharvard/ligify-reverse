from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pandas as pd
import re
import ast
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent / "src"))

from chemical import get_smiles_from_chebi

def tanimoto_calc(smi1, smi2):
    try:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s
    except:
        return 0

# try newer scorings
def tanimoto_calc_new(smi1, smi2):
    try:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        mfpgen = AllChem.GetMorganGenerator(radius=3, fpSize = 2048)
        fp1 = mfpgen.GetFingerprint(mol1)
        fp2 = mfpgen.GetFingerprint(mol2)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s
    except:
        return 0


def getChemicalsFromLogFile(acc):
    folder_path = '/Users/jiojeong/Documents/ligify-reverse-1/benchmarking/log_files/log remote blast/tanimoto/'
    log_file_path = folder_path + str(acc) + ".log"
    with open(log_file_path, 'r') as file:
        content = file.read()
    
    # Find the starting and ending indices of the chemicals section
    start_marker = "These are all the chemicals ranked:"
    end_marker = "________________________________"
    start_index = content.find(start_marker) + len(start_marker)
    end_index = content.find(end_marker, start_index)
    
    # Extract the section that contains the chemical entries
    chemicals_section = content[start_index:end_index].strip()
    
    # Split the section into individual lines, each representing a chemical dictionary
    chemical_entries = chemicals_section.split('\n')
    
    # Convert each string representation of a dictionary to an actual dictionary
    chemicals = [ast.literal_eval(entry) for entry in chemical_entries if entry.strip()]
    
    return chemicals

def add_smiles_to_dict(chem_list):
    for chemical in chem_list:
        if (type(chemical) is list):
            smiles = get_smiles_from_chebi(chemical[2])
            chemical.append(smiles)
        else:
            smiles = get_smiles_from_chebi(chemical["ChEBI ID"])
            chemical['SMILES'] = smiles
    return chem_list
        

def tanimoto_log(true_ligand_smiles, chem_list):
    if (chem_list == []):
        return []
    for chemical in chem_list:
        compare = ""
        if (type(chemical) is list):
            compare = chemical[-1]
        else:
            compare = chemical['SMILES']
        #tanimoto = tanimoto_calc(true_ligand_smiles, compare)
        tanimoto_new = tanimoto_calc_new(true_ligand_smiles, compare)
        if (type(chemical) is list):
            #chemical.append(tanimoto)
            chemical.append(tanimoto_new)
        else:
            #chemical["Tanimoto with true ligand"] = tanimoto
            chemical["Tanimoto with true ligand"] = tanimoto_new
    
    return chem_list
        
      
    
if __name__ == "__main__":    
    excel_file_path = "/Users/jiojeong/Documents/ligify-reverse-1/benchmarking/Benchmarking_Dataset.xlsx"
    df = pd.read_excel(excel_file_path)

    # Assuming column headers are present and you want to use columns A and D
    # If your columns don't have headers, you can specify them using 'header=None' parameter in read_excel() and provide 'names' parameter to assign column names.

    # Create the dictionary
    excel_dict = dict(zip(df['Biosensor RefSeq'], df['SMILES']))
    
    queries = ["NP_414606.1", "WP_004925500.1", "WP_011336736.1", "NP_414847.3", "NP_418026.1", "NP_415533.1", "WP_011014162.1", "NP_862536.1", "NP_418342.2", "WP_012368846.1", "WP_010974118.1", "NP_418209.1", "NP_414879.3", "WP_000174305.1", "NP_391277.1", "AAC37157.1", "BAE47075.1", "NP_746731.1", "WP_002965779.1", "WP_009968057.1", "WP_011594778.1", "NP_745052.1", "WP_003227022.1", "WP_010813722.1", "WP_001278727.1", "NP_415039.1", "WP_010811303.1", "WP_003114242.1", "WP_013056567.1", "NP_414655.1", "WP_011731512.1",
    'WP_003183369.1', 'WP_010811734.1', 'WP_011614884.1', 'WP_011615200.1', 'WP_010810611.1',
    'WP_010812428.1', 'NP_421195.1', 'WP_010809162.1', 'BAA86295.1', 'AAK15050.1', 'WP_011614439.1',
    'WP_041688483.1', 'WP_010951546.1', 'WP_013233032.1', 'WP_010813655.1', 'WP_011616377.1',
    'WP_001137892.1', 'WP_000929443.1', 'WP_000234823.1', 'WP_001300658.1', 'WP_001300658.1',
    'WP_035269502.1', 'WP_009944749.1', 'ACS29497.1', 'AAB62296.1', 'WP_000113282.1', 'BAA03510.1',
    'NP_266817.1', 'NP_059715.1', 'NP_417391.1', 'NP_396546.3', 'WP_004926797.1', 'NP_742225.1',
    'NP_721604.1', 'NP_542858.1', 'NP_635750.1', 'ACY33523.1', 'BAK65962.1', 'NP_630777.1',
    'NP_416175.1', 'WP_011291385.1', 'WP_034117836.1', 'WP_002397724.1', 'ABP48117.1', 'BAH70273.1',
    'WP_014076817.1', 'WP_142664765.1', 'AAK38101.1', 'WP_011028828.1', 'NP_418810.1', 'A0A0D5A3S5',
    'NP_059701.1', 'WP_003514478.1', 'APY21445.1', 'WP_001145439.1', 'WP_011156222.1', 'WP_004926705.1',
    'WP_051824537.1', 'WP_011030045.1', 'WP_015475612.1', 'WP_079651080.1', 'BAK67179.1', 'AAA25771.1',
    'AAW51730.1', 'AUW46298.1', 'WP_012273540.1', 'WP_011004209.1', 'AAY86547.1', 'CAY46636.1'
    ]
    

    # Print the dictionary to verify
    for acc in queries:
        if (acc == "nan"):
            continue
        
        print("Running for: " + acc)
        log_file_name = acc + '.log'
        with open(log_file_name, 'w') as log_file:
            sys.stdout = log_file
        
            print("Query: " + acc)
            chemicals = getChemicalsFromLogFile(acc)
            print("Chemicals: ")
            for i in chemicals:
                print(i)
            true_ligand_smiles = excel_dict[acc]
            
            print("true ligand smiles: ")
            print(true_ligand_smiles)
            
            chemicals = add_smiles_to_dict(chemicals)
            modified_chem_list = tanimoto_log(true_ligand_smiles, chemicals)
            print("Modified chemicals list: ")
            for i in modified_chem_list:
                print(i)
            
            sys.stdout = sys.__stdout__