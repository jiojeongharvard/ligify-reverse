import subprocess
import requests
import json
import pandas as pd
import streamlit as st
from accID2operon import acc2operon
from blast import blast
import sys

import xml.etree.ElementTree as ET
import re
import datetime
import time

# filter out the dataframe of homologs based on given cutoffs
def filterHomologs(homologs: pd.DataFrame, ident_cutoff:int , cov_cutoff: int):
    filtered_df = homologs[(homologs['Identity'] > ident_cutoff) & (homologs['Coverage'] > cov_cutoff)]
    filtered_df = filtered_df[~((filtered_df['Identity'] == 100) & (filtered_df['Coverage'] == 100))]
    return filtered_df

# use the acc2operon to get the operons for the input biosensor and its homologs
def getOperon(acc: str, homolog_df: pd.DataFrame):
    input_reg_operon = acc2operon(acc)
    homolog_df['NCBI Id'] = homolog_df['Uniprot Id'].apply(lambda x: get_ncbi_id(x))
    homolog_df['EMBL Id'] = homolog_df['Uniprot Id'].apply(lambda x: convert_uniprot_to_embl(x))
    homolog_df['Operon'] = homolog_df['NCBI Id'].apply(lambda x: acc2operon(x))
    #homolog_df['Operon'] = homolog_df['EMBL Id'].apply(lambda x: acc2operon(x))
    
    homolog_df = homolog_df[homolog_df['Operon'] != 'EMPTY']
    
    # homolog_df['duplicate test'] = homolog_df['Operon'].apply(lambda x: x['operon'][0]['alias'])
    # homolog_df.drop_duplicates(subset='duplicate test', keep='first', inplace=True)
    
    homolog_df.drop_duplicates(subset='Operon', keep='first', inplace=True)
    
    # filter_condition = homolog_df['Operon'].apply(lambda x: x['operon'][0]['alias']) == homolog_df['Operon'].apply(lambda x: x['operon'][0]['alias']).iloc[0]
    # homolog_df = homolog_df[~filter_condition].reset_index(drop=True)
    #homolog_df.drop_duplicates(subset='duplicate check', keep='first', inplace=True)
    # homolog_df['Operon'] = homolog_df['Operon'].apply(lambda x: list({d['alias']: d for d in x['operon']}.values()))
    
    return input_reg_operon, homolog_df

# take out / delete all the regulators from a given operon
def filterRegFromOperon(acc: str, input_reg_operon):
    index = 0
    index_of_reg = None
    for gene in input_reg_operon["operon"]:
        gene['index'] = index
        if gene['accession'] == acc:
            index_of_reg = index
        index = index + 1
    
    if index_of_reg != None:
        for gene in input_reg_operon["operon"]:
            if gene['accession'] != acc:
                gene['distance'] = abs(gene['index'] - index_of_reg)
    
    regulator = re.compile(r"regulator|repressor|activator")
    
    to_replace = []
    for gene in input_reg_operon["operon"]:
        if "description" in gene.keys():
            if not regulator.search(gene["description"]):
                to_replace.append(gene)
                
    input_reg_operon["operon"] = to_replace
     
    return input_reg_operon

# take out / delete all the regulators from a given dataframe of operons
def filterRegFromOperondf(homolog_operon_df: pd.DataFrame):
    for index, row in homolog_operon_df.iterrows():
        if row['Operon'] == 'EMPTY':
            continue
        else:
            acc = row["NCBI Id"]
            input_operon = row["Operon"]
            filterRegFromOperon(acc, input_operon)
    
    return homolog_operon_df

# this is the exact same function as the protein2chemicals function in regular ligify
def getChemicalsFromAcc(acc):    
    url = "https://rest.uniprot.org/uniprotkb/search?query="+acc+"&format=json"

    response = requests.get(url)
    if response.ok:
        protein = json.loads(response.text)


        if len(protein['results']) > 0:
            if "comments" in protein["results"][0]:
                protein = protein["results"][0]["comments"]

                protein_data = {}

                    # look for description
                Function = [i["texts"][0]["value"] for i in protein if i["commentType"] == "FUNCTION"]
                if len(Function) != 0:
                    protein_data["function"] = Function[0]
        

                    # look for catalytic activity
                CATALYSIS = [i["reaction"]["name"] for i in protein if i["commentType"] == "CATALYTIC ACTIVITY"]
                if len(CATALYSIS) != 0:
                    protein_data["catalysis"] = CATALYSIS[0]
                    

                    # look for catalytic activity
                try:
                    RXN = [i["reaction"]["reactionCrossReferences"] for i in protein if i["commentType"] == "CATALYTIC ACTIVITY"]
                    if len(RXN) != 0:
                        LIGANDS = [i["id"] for i in RXN[0] if i["database"] == "ChEBI"]
                        if len(LIGANDS) != 0:
                            protein_data["ligands"] = LIGANDS
                except:
                    pass

                    # look for induction
                INDUCTION = [i["texts"][0]["value"] for i in protein if i["commentType"] == "INDUCTION"]
                if len(INDUCTION) != 0:
                    protein_data["induction"] = INDUCTION[0]

                            # look for induction
                PATHWAY = [i["texts"][0]["value"] for i in protein if i["commentType"] == "PATHWAY"]
                if len(PATHWAY) != 0:
                    protein_data["pathway"] = PATHWAY[0]

                return(protein_data)
    else:
        response.raise_for_status()
    
# convert chebi id to the chemical name
def chebi_id_to_name(chebi_id):
    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
    # Send an HTTP request to the URL
    response = requests.get(url)

    # Extract the HTML content
    html_content = response.text

    # Define a regular expression pattern to match the ChEBI Name within the <title> tag
    pattern = r'<title>\s*(.*?)\s*\(\s*CHEBI:\d+\s*\)\s*</title>'

    # Search for the pattern in the HTML content
    match = re.search(pattern, html_content)

    # Check if a match is found
    if match:
        # Extract the ChEBI Name
        chebi_name = match.group(1).strip()
        return chebi_name
    else:
        print("ChEBI Name not found.")
        return None
    
    
# iterate through the proteins in the operon to look for chemicals in the reaction database. 
# add the list of chemicals to a new field named "chemicals" for each protein
# thus each "operon" will have a list of "genes/proteins" that each have a list of "chemicals"       
def addChemicalsToOperon(input_reg_operon):
    not_ligands = ["hydron", "water", "carbon dioxide", 'NADP(3-)', 'NADPH(4-)', 'NADH(2-)', 'NAD(1-)', 'hydrogen peroxide', 'dioxygen', "coenzyme A(4-)",
                    "acetyl-CoA(4-)","acyl-CoA(4-)","hydrogenphosphate","ATP(4-)","ADP(3-)","diphosphate(3-)","DNA 5'-phosphate polyanion","hydrogen acceptor",
                    "hydrogen donor","FMNH(2)(2-)","FMN(3-)","potassium(1+)","iron(2+)","iron(3+)","copper(1+)","copper(2+)","ammonium","aldehyde","carboxylic acid anion",
                    "superoxide","AMP 3'-end(1-) residue","thiol group","H group","lipid II(3-)"]
    for gene in input_reg_operon["operon"]:
        ligand_list = []
        protein_data = getChemicalsFromAcc(gene["accession"])
        if isinstance(protein_data, dict):
            if "ligands" in protein_data.keys():
                ligands = protein_data["ligands"]
                for ligand in ligands:
                    ligand_name = chebi_id_to_name(ligand)
                    if ligand_name != None and ligand_name not in not_ligands:
                        ligand_dict = {"name": ligand_name, "ChEBI ID": ligand}
                        ligand_list.append(ligand_dict)
                        
        if ligand_list:
            gene['chemicals'] = ligand_list
        
        
    
    return input_reg_operon
    
# adding chemicals to a DATAFRAME of multiple operons     
def addChemicalsToOperondf(homolog_operon_df: pd.DataFrame):
    for index, row in homolog_operon_df.iterrows():
        if row['Operon'] == 'EMPTY':
            continue
        else:
            input_operon = row["Operon"]
            addChemicalsToOperon(input_operon)
    
    return homolog_operon_df
    
# makes a dictionary where the keys are the chemical ids.
# the value is a list of all the times the chemical appears in the operons   
def getAllChemicalsInOperons(input_reg_operon, homolog_operon_df):
    all_chemicals = {}
    chemicals_in_input_reg_operon = []
    total_gene_num = 0
    
    for gene in input_reg_operon["operon"]:
        if "chemicals" in gene.keys():
            total_gene_num += 1
            for chemical in gene["chemicals"]:
                id = chemical["ChEBI ID"]
                
                info = {}
                if "distance" in gene.keys():
                    info["distance_from_reg"] = gene["distance"]
                info["identity"] = 100
                info["coverage"] = 100
                
                all_chemicals.setdefault(id, []).append(info)
                chemicals_in_input_reg_operon.append(id)
    
    if homolog_operon_df is None:
        return all_chemicals, set(chemicals_in_input_reg_operon), total_gene_num
    
    for index, row in homolog_operon_df.iterrows():
        op = row['Operon']
        if op != 'EMPTY':
            for gene in op['operon']:
                if "chemicals" in gene.keys():
                    total_gene_num += 1
                    for chemical in gene["chemicals"]:
                        id = chemical["ChEBI ID"]
                        
                        info = {}
                        info = {}
                        if "distance" in gene.keys():
                            info["distance_from_reg"] = gene["distance"]
                        info["identity"] = row["Identity"]
                        info["coverage"] = row["Coverage"]
                        
                        all_chemicals.setdefault(id, []).append(info)        
            
    return all_chemicals, set(chemicals_in_input_reg_operon), total_gene_num

# compute the score of the chemical based on identity and coverage of homolog and distance between GOI and regulator
# weights can be adjusted
def rankChemicals(all_chemicals, chemicals_in_input_reg_operon, cutoff, num_genes, weights):
    chemical_scores = []
    
    distance_weight = weights["distance"]
    identity_weight = weights["identity"]
    coverage_weight = weights["coverage"]
    num_occurrence_weight = weights["frequency"]
    query_operon_max_penalty = weights["max input operon penalty"]
    
    for chemical in all_chemicals:
        distance_counter = 0
        distance_sum = 0
        identity_sum = 0
        coverage_sum = 0
        # homolog similarities > how frequent > distance
        
        occurrences = all_chemicals[chemical]
        num_occurrence = len(occurrences)
        
        for occurrence in occurrences:
            if "distance_from_reg" in occurrence:
                distance_counter += 1
                distance_sum += occurrence['distance_from_reg']
            identity_sum += occurrence['identity']
            coverage_sum += occurrence['coverage']
            
        identity_avg = identity_sum / num_occurrence    
        coverage_avg = coverage_sum / num_occurrence
        
        min_occurrence = 0
        max_occurrence = num_genes
        normalized_occurrence = (num_occurrence - min_occurrence) / (max_occurrence - min_occurrence)
        # print("chemical id: " + str(chemical))
        # print("identity avg: " + str(identity_avg) + "; coverage avg: " + str(coverage_avg) + "; normalized occ: " + str(normalized_occurrence))
        
        not_present_in_query_operon_penalty = 0
        if chemical not in chemicals_in_input_reg_operon:
            not_present_in_query_operon_penalty = query_operon_max_penalty

        
        # print("not in query operon penalty: " + str(not_present_in_query_operon_penalty))   
        
        sub_scores = {}
        if distance_counter > 0:
            distance_avg = distance_sum / distance_counter
            # print("distance avg: " + str(distance_avg))
            total_score = (1000 - not_present_in_query_operon_penalty - (distance_avg * distance_weight) - (100 - identity_avg) * identity_weight - (100 - coverage_avg) * coverage_weight - (1 - normalized_occurrence) * num_occurrence_weight) / 10
            sub_scores["query operon penalty"] = not_present_in_query_operon_penalty
            sub_scores["distance average"] = distance_avg
            sub_scores["identity average"] = identity_avg
            sub_scores["coverage average"] = coverage_avg
            sub_scores["frequency"] = normalized_occurrence
        else:
            total_score =  (1000 - not_present_in_query_operon_penalty - (distance_weight * 5) - (100 - identity_avg) * identity_weight - (100 - coverage_avg) * coverage_weight - (1 - normalized_occurrence) * num_occurrence_weight) / 10
            sub_scores["query operon penalty"] = not_present_in_query_operon_penalty
            sub_scores["identity average"] = identity_avg
            sub_scores["coverage average"] = coverage_avg
            sub_scores["frequency"] = normalized_occurrence
        
        if (total_score > cutoff):
            chemical_scores.append([chemical, total_score, chemical, sub_scores])
            
    sorted_chemical_scores = sorted(chemical_scores, key=lambda x: x[1], reverse = True)
    
    for l in sorted_chemical_scores:
        l[0] = chebi_id_to_name(l[0])
        
    
    num_occurrence_weight = 150
    return sorted_chemical_scores           
            
      
def convert_uniprot_to_embl(uniprot_id):
    base_url = "https://www.uniprot.org/uniprot/"
    url = base_url + uniprot_id + ".xml"
    response = requests.get(url)
    if response.status_code == 200:
        # Parse XML response to get RefSeq ID
        # xml_data = response.text
        # print(xml_data)
        # root = ET.fromstring(xml_data)
        
        # embl_id = None
        # for dbReference in root.findall(".//{http://uniprot.org/uniprot}dbReference"):
        #     if dbReference.attrib['type'] == 'EMBL':
        #         embl_id = dbReference.attrib['id']
        #         break
        # return embl_id
        
        # Parse XML
        xml_data = response.text
        root = ET.fromstring(xml_data)

        # Define the namespace
        ns = {'uni': 'http://uniprot.org/uniprot'}

        # Find the dbReference element with type="EMBL"
        embl_element = root.find(".//uni:dbReference[@type='EMBL']", ns)

        # Extract the value of the property element where type="protein sequence ID"
        if embl_element is not None:
            protein_sequence_id = embl_element.find("./uni:property[@type='protein sequence ID']", ns)
            if protein_sequence_id is not None:
                embl_id = protein_sequence_id.attrib.get('value')
                # print("EMBL ID:", embl_id)
                return embl_id
            else:
                print("Protein sequence ID not found.")
        else:
            print("EMBL element not found.")

    else:
        print("Error: Unable to fetch data from UniProt")
        return None


# converting from refseq -> uniprot id
def get_uniprot_id(refseq_accession):
    url = "https://rest.uniprot.org/uniprotkb/search?query="+refseq_accession+"&format=json"
    response = requests.get(url)

    if response.status_code == 200:
        data = json.loads(response.text)["results"][0]
        acc = data["primaryAccession"]
        return acc

    return None
    
# converting from uniprot id -> ncbi id
def get_ncbi_id(uniprot_id):
    base_url = "https://www.uniprot.org/uniprot/"
    url = base_url + uniprot_id + ".xml"
    response = requests.get(url)
    if response.status_code == 200:
        # Parse XML response to get RefSeq ID
        xml_data = response.text
        root = ET.fromstring(xml_data)
        refseq_id = None
        for dbReference in root.findall(".//{http://uniprot.org/uniprot}dbReference"):
            if dbReference.attrib['type'] == 'RefSeq':
                refseq_id = dbReference.attrib['id']
                break
        return refseq_id
    else:
        print("Error: Unable to fetch data from UniProt")
        return None

def get_enzyme_description(input_operon):
    all_enzymes_list = []
    for gene in input_operon["operon"]:
        if "description" in gene.keys():
            all_enzymes_list.append(gene["description"])
            
    return all_enzymes_list

def get_enzyme_description_df(homolog_operon_df: pd.DataFrame):
    all_enzymes_list = []
    for index, row in homolog_operon_df.iterrows():
        if row['Operon'] == 'EMPTY':
            continue
        else:
            input_operon = row["Operon"]
            temp = get_enzyme_description(input_operon)
            all_enzymes_list.extend(temp)
    
    return all_enzymes_list
  
    
def log_function(query, excel_dict):
    print("Running for " + query)
    start_time = datetime.datetime.now()
    # Create or open the log file with the name as the input string
    log_file_name = query + '.log'
    
    with open(log_file_name, 'w') as log_file:
        # Redirect standard output to the log file
        sys.stdout = log_file
        params = {"ident_cutoff": 70, "cov_cutoff": 90}
        df = blast(query, "RefSeq", params, 20)
        filtered_df = filterHomologs(df, 70, 90)
        # too lenient 90 percent coverage 60-70 identity      
        print("Query: " + query)
        inputs_operon, homolog_operons = getOperon(query, filtered_df)
        print("This is how the query's operon looks like:")
        print(inputs_operon)
        print("________________________________")
        
        print("This is how the dataframe of homologs looks like:")
        print(homolog_operons)
        print("________________________________")
        
        # print("This is how the homologs' operons look like:")
        # for index, row in homolog_operons.iterrows():
        #     op = row['Operon']
        #     if op != 'EMPTY':
        #         print(op)
        # print("________________________________")
        
        filterRegFromOperon(query, inputs_operon)
        filterRegFromOperondf(homolog_operons)
        
        print("This contains the description of all enzymes (excluding regulators) in the operons")
        temp1 = get_enzyme_description(inputs_operon)
        temp2 = get_enzyme_description_df(homolog_operons)
        temp1.extend(temp2)
        print(temp1)
        print("________________________________") 
        
        addChemicalsToOperon(inputs_operon)
        addChemicalsToOperondf(homolog_operons)
        
        print("This is how the chemicals are added to the genes in the operon:")
        for gene in inputs_operon["operon"]:
            if "chemicals" in gene.keys():
                print(gene["chemicals"])
        print("________________________________")       
        
        print("This is how the chemicals are added to the genes in the homolog operons:")        
        for index, row in homolog_operons.iterrows():
            op = row['Operon']
            if op != 'EMPTY':
                for gene in op['operon']:
                    if "chemicals" in gene.keys():
                        print(gene['chemicals'])
        print("________________________________")   
                      
        
        print("These are all the chemicals in the operons:")
        chem, chemicals_in_input_reg_operon, num_genes = getAllChemicalsInOperons(inputs_operon, homolog_operons)
        print(chem)
        print("________________________________")      
        
        print("The total number of enzymes that have chemicals found in the query operon + all homolog operons is: " + str(num_genes))  
        print("________________________________")  
        
        print("These are all the chemicals ranked:")
        weights = {}
        weights["distance"] = 5
        weights["identity"] = 11
        weights["coverage"] = 20
        weights["frequency"] = 170
        weights["max input operon penalty"] = 250
        output_chemicals, used_weights = rankChemicals(chem, chemicals_in_input_reg_operon, 0, num_genes, weights)
        for i in output_chemicals:
            print(i)
        
        print("________________________________")  
        print("This is the correct chemical for the query: ")
        print(excel_dict[query])
           
        print("________________________________")     
        
        print("These are weights used for ranking:")
        print(weights)
    
        
        print("________________________________")      
        
        end_time = datetime.datetime.now()
        runtime = end_time - start_time
        print("The Runtime of query: " + str(runtime))

        # Your function logic that prints statements
        
    # Restore standard output
    sys.stdout = sys.__stdout__
      
if __name__ == "__main__":
    #all queries that work with ligify
    all_queries = ["NP_414606.1", "WP_004925500.1", "WP_011336736.1", "NP_414847.3", "NP_418026.1", "NP_415533.1", "WP_011014162.1", "NP_862536.1", "NP_418342.2", "WP_012368846.1", "WP_010974118.1", "NP_418209.1", "NP_414879.3", "WP_000174305.1", "NP_391277.1", "AAC37157.1", "BAE47075.1", "NP_746731.1", "WP_002965779.1", "WP_009968057.1", "WP_011594778.1", "NP_745052.1", "WP_003227022.1", "WP_010813722.1", "WP_001278727.1", "NP_415039.1", "WP_010811303.1", "WP_003114242.1", "WP_013056567.1", "NP_414655.1", "WP_011731512.1"]
    
    #queries that work in ligify + have no homolog hits
    no_homolog = ["NP_414606.1", "NP_414655.1", "NP_414879.3", "NP_415039.1", "NP_415533.1", "NP_418026.1", "NP_418209.1", "WP_000174305.1", "WP_001278727.1", "WP_002965779.1", "WP_003114242.1", "WP_013056567.1"]
    
    #queries that work in ligify + have at least a few homolog hits
    queries = [elem for elem in all_queries if elem not in no_homolog]
    
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
    
    # Read the Excel file
    excel_file_path = '/Users/jiojeong/Desktop/Benchmarking_Dataset.xlsx'
    df = pd.read_excel(excel_file_path)

    # Assuming column headers are present and you want to use columns A and D
    # If your columns don't have headers, you can specify them using 'header=None' parameter in read_excel() and provide 'names' parameter to assign column names.

    # Create the dictionary
    excel_dict = dict(zip(df['Biosensor RefSeq'], df['Chemical name']))

    # Print the dictionary to verify
    print(excel_dict)
    
    for query in queries:
        log_function(query, excel_dict)
    
    # print(acc2operon("WP_016167756.1"))
    # print("__")
    # print(acc2operon("AUX87632.1"))
    # print("__")
    # print(acc2operon("TCB33756.1"))
    
    # oper = acc2operon("WP_016167756.1")
    # print(oper)
    # print(oper['operon'][0]['alias'])
    
    # data = {
    # 'Name': ['A', 'B', 'C', 'D', 'E', 'F'], 
    # 'Operon': [{'operon':[{'alias': "two"}, {'alias': "wrong"}], 'thing': "okay"}, 
	# 	{'operon':[{'alias': "one"}, {'alias': "wrong"}], 'thing': "okay"}, 
    #            {'operon':[{'alias': "one"}, {'alias': "wrong"}], 'thing': "okay"}, 
    #            {'operon':[{'alias': "three"}, {'alias': "wrong"}], 'thing': "okay"}, 
	# 	{'operon':[{'alias': "three"}, {'alias': "wrong"}], 'thing': "okay"},
	# 	{'operon':[{'alias': "two"}, {'alias': "wrong"}], 'thing': "okay"}] 
    # }
    # # Creating DataFrame
    # dff = pd.DataFrame(data)
    # print(dff)
    # # filter_condition = dff['Operon'].apply(lambda x: x['operon'][0]['alias']) == dff['Operon'].apply(lambda x: x['operon'][0]['alias']).iloc[0]
    # # dff = dff[~filter_condition].reset_index(drop=True)
    # dff['test'] = dff['Operon'].apply(lambda x: x['operon'][0]['alias'])
    # dff.drop_duplicates(subset='test', keep='first', inplace=True)

    # print(dff)
    
