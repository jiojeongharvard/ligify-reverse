import streamlit as st

import pandas as pd

import sys
from pathlib import Path

# Append the parent directory of "src" to the system path
sys.path.append(str(Path(__file__).parent.parent / "src"))

from operon import log_function

if __name__ == "__main__":
    #all queries that work with ligify
    all_queries = ["NP_414606.1", "WP_004925500.1", "WP_011336736.1", "NP_414847.3", "NP_418026.1", "NP_415533.1", "WP_011014162.1", "NP_862536.1", "NP_418342.2", "WP_012368846.1", "WP_010974118.1", "NP_418209.1", "NP_414879.3", "WP_000174305.1", "NP_391277.1", "AAC37157.1", "BAE47075.1", "NP_746731.1", "WP_002965779.1", "WP_009968057.1", "WP_011594778.1", "NP_745052.1", "WP_003227022.1", "WP_010813722.1", "WP_001278727.1", "NP_415039.1", "WP_010811303.1", "WP_003114242.1", "WP_013056567.1", "NP_414655.1", "WP_011731512.1"]
    
    #queries that work in ligify + have no homolog hits
    no_homologs = ["NP_414606.1", "NP_414655.1", "NP_414879.3", "NP_415039.1", "NP_415533.1", "NP_418026.1", "NP_418209.1", "WP_000174305.1", "WP_001278727.1", "WP_002965779.1", "WP_003114242.1", "WP_013056567.1"]
    
    #queries that work in ligify + have at least a few homolog hits
    yes_homologs = [elem for elem in all_queries if elem not in no_homologs]
    
    # these are all the queries that don't work with ligify
    failure_mode_queries = [
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
    
    queries =  ["WP_011731512.1"]
    
    # Read the Excel file
    excel_file_path = Path(__file__).parent / "Benchmarking_Dataset.xlsx"
    df = pd.read_excel(excel_file_path)

    # Assuming column headers are present and you want to use columns A and D
    # If your columns don't have headers, you can specify them using 'header=None' parameter in read_excel() and provide 'names' parameter to assign column names.

    # Create the dictionary of correct biosensor-chemical pairs
    excel_dict = dict(zip(df['Biosensor RefSeq'], df['Chemical name']))

    # Print the dictionary to verify
    print(excel_dict)
    
    # Run the log function for each query
    for query in queries:
        log_function(query, excel_dict)