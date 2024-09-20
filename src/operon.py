import requests
import json
import pandas as pd
import streamlit as st
from accID2operon import acc2operon
from blast import blast, blast_remote
import sys

import xml.etree.ElementTree as ET
import re
import datetime
from collections import Counter

# filter out the dataframe of homologs based on given cutoffs
def filterHomologs(homologs: pd.DataFrame, ident_cutoff:int , cov_cutoff: int):
    if (homologs is None):
        return None
    filtered_df = homologs[(homologs['Identity'] > ident_cutoff) & (homologs['Coverage'] > cov_cutoff)]
    filtered_df = filtered_df[~((filtered_df['Identity'] == 100) & (filtered_df['Coverage'] == 100))]
    return filtered_df

# use the acc2operon functino to get the operons for the input biosensor and its homologs
def getOperon(acc: str, homolog_df: pd.DataFrame):
    input_reg_operon = acc2operon(acc)
    
    if (homolog_df is None):
        return input_reg_operon, None
    
    if 'Uniprot Id' in homolog_df.columns:
        homolog_df['NCBI Id'] = homolog_df['Uniprot Id'].apply(lambda x: get_ncbi_id(x))
        homolog_df['EMBL Id'] = homolog_df['Uniprot Id'].apply(lambda x: convert_uniprot_to_embl(x))
    
    if 'NCBI Id' in homolog_df.columns:
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
    if (input_reg_operon == "EMPTY"):
        return "EMPTY"
    
    regulator = re.compile(r"regulator|repressor|activator|regulatory")
    
    index = 0
    index_of_reg = input_reg_operon["enzyme_index"]
    
    for gene in input_reg_operon["operon"]:
        gene['index'] = index
        index = index + 1
    
    if index_of_reg != None:
        for gene in input_reg_operon["operon"]:
            gene['distance'] = abs(gene['index'] - index_of_reg)
    
    to_replace = []
    for gene in input_reg_operon["operon"]:
        if "description" in gene.keys():
            if not regulator.search(gene["description"]) and gene['accession'] != acc:
                to_replace.append(gene)
                
    input_reg_operon["operon"] = to_replace
     
    return input_reg_operon

# take out / delete all the regulators from a given dataframe of operons
def filterRegFromOperondf(homolog_operon_df: pd.DataFrame):
    if (homolog_operon_df is None):
        return None
    
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
                    protein_data["numReactions"] = len(RXN)
                except:
                    protein_data["numReactions"] = 0
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
    if (input_reg_operon == 'EMPTY'):
        return "EMPTY"
    
    not_ligands = ["hydron", "water", "carbon dioxide", 'NADP(3-)', 'NADPH(4-)', 'NADH(2-)', 'NAD(1-)', 'hydrogen peroxide', 'dioxygen', "coenzyme A(4-)",
                    "acetyl-CoA(4-)","acyl-CoA(4-)","hydrogenphosphate","ATP(4-)","ADP(3-)", "AMP(2-)", "diphosphate(3-)","DNA 5'-phosphate polyanion","hydrogen acceptor",
                    "hydrogen donor","FMNH(2)(2-)","FMN(3-)","potassium(1+)","iron(2+)","iron(3+)","copper(1+)","copper(2+)","ammonium","aldehyde","carboxylic acid anion",
                    "superoxide","AMP 3'-end(1-) residue","thiol group","H group","lipid II(3-)", 'GTP(4-)', 'GDP(3-)', 'GMP(2-)' ,'FAD(3-)', 'FADH(2)(2-)', ' sodium(1+)', 
                    'flavin(1-)', '1,5-dihydroflavin']
    
    reaction_count = 0
     
    for gene in input_reg_operon["operon"]:
        ligand_list = []
        protein_data = getChemicalsFromAcc(gene["accession"])
        if isinstance(protein_data, dict):
            reaction_count += protein_data["numReactions"]
            if "ligands" in protein_data.keys():
                ligands = protein_data["ligands"]
                for ligand in ligands:
                    ligand_name = chebi_id_to_name(ligand)
                    if ligand_name != None and ligand_name not in not_ligands:
                        ligand_dict = {"name": ligand_name, "ChEBI ID": ligand}
                        ligand_list.append(ligand_dict)
                        
        if ligand_list:
            gene['chemicals'] = ligand_list
        
        
    input_reg_operon["totalRxnCount"] = reaction_count
    return input_reg_operon
    
# adding chemicals to a DATAFRAME of multiple operons     
def addChemicalsToOperondf(homolog_operon_df: pd.DataFrame):
    if (homolog_operon_df is None):
        return None
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
    total_enz_num = 0
    enz_w_lig_num = 0
    total_rxn_num = 0
    
    if (input_reg_operon != "EMPTY"):
        total_rxn_num += input_reg_operon["totalRxnCount"]
        for gene in input_reg_operon["operon"]:
            total_enz_num += 1
            if "chemicals" in gene.keys():
                enz_w_lig_num += 1
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
        return all_chemicals, set(chemicals_in_input_reg_operon), total_enz_num, enz_w_lig_num, total_rxn_num
    
    for index, row in homolog_operon_df.iterrows():
        op = row['Operon']
        if op != 'EMPTY':
            total_rxn_num += op["totalRxnCount"]
            for gene in op['operon']:
                total_enz_num += 1
                if "chemicals" in gene.keys():
                    enz_w_lig_num += 1
                    for chemical in gene["chemicals"]:
                        id = chemical["ChEBI ID"]
                        
                        info = {}
                        info = {}
                        if "distance" in gene.keys():
                            info["distance_from_reg"] = gene["distance"]
                        info["identity"] = row["Identity"]
                        info["coverage"] = row["Coverage"]
                        
                        all_chemicals.setdefault(id, []).append(info)        
            
    return all_chemicals, set(chemicals_in_input_reg_operon), total_enz_num, enz_w_lig_num, total_rxn_num

# compute the score of the chemical based on identity and coverage of homolog and distance between GOI and regulator
# weights can be adjusted
def rankChemicals(all_chemicals, chemicals_in_input_reg_operon, cutoff, enz_w_lig_num, weights):
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
        max_occurrence = enz_w_lig_num
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
            temp_dict = {"Chemical Name": chemical, "Score": total_score, "ChEBI ID": chemical, "Subscore": sub_scores}
            chemical_scores.append(temp_dict)
            
    sorted_chemical_scores_list_of_dict = sorted(chemical_scores, key=lambda x: x["Score"], reverse = True)

    for l in sorted_chemical_scores_list_of_dict:
        l["Chemical Name"] = chebi_id_to_name(l["ChEBI ID"])
        
    return sorted_chemical_scores_list_of_dict           
            
      
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

    if response.ok:
        if (json.loads(response.text)["results"] == []):
            print("FATAL: Bad uniprot API request "+ str(response.status_code))
            st.error("RefSeq cannot be found in UniprotKB")
            return ""
        data = json.loads(response.text)["results"][0]
        acc = data["primaryAccession"]
        return acc
    else:
        print("FATAL: Bad uniprot API request "+ str(response.status_code))
        st.error("RefSeq cannot be found in UniprotKB")
        return ""
    
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

# get the description of the enzymes in the operon
def get_enzyme_description(input_operon):
    all_enzymes_list = []
    if (input_operon == "EMPTY"):
        return all_enzymes_list
    
    for gene in input_operon["operon"]:
        if "description" in gene.keys():
            all_enzymes_list.append(gene["description"])
            
    return all_enzymes_list

# get the description of the enzymes in a dataframe of operons
def get_enzyme_description_df(homolog_operon_df: pd.DataFrame):
    all_enzymes_list = []
    if (homolog_operon_df is None):
        return all_enzymes_list
    for index, row in homolog_operon_df.iterrows():
        if row['Operon'] == 'EMPTY':
            continue
        else:
            input_operon = row["Operon"]
            temp = get_enzyme_description(input_operon)
            all_enzymes_list.extend(temp)
    
    return all_enzymes_list

# used for analyzing enzyme descriptions
def list_of_phrases_to_frequencies(list_of_phrases):
    phrase_counts = Counter(list_of_phrases)
    sorted_phrases = dict(sorted(phrase_counts.items(), key=lambda item: item[1], reverse=True))
    return sorted_phrases


def list_of_phrases_to_word_to_frequencies(list_of_phrases):
    list_of_words = [word for phrase in list_of_phrases for word in phrase.split()]
    word_counts = Counter(list_of_words)
    sorted_words = dict(sorted(word_counts.items(), key=lambda item: item[1], reverse=True))
    return sorted_words

# log function to look at all the relevant data    
def log_function(query, excel_dict):
    print("Running for " + query)
    start_time = datetime.datetime.now()
    # Create or open the log file with the name as the input string
    log_file_name = query + '.log'
    
    with open(log_file_name, 'w') as log_file:
        # Redirect standard output to the log file
        sys.stdout = log_file
        print("Query: " + query)
        params = {"ident_cutoff": 70, "cov_cutoff": 90}
        df = blast_remote(query, "RefSeq", params, 20)
        filtered_df = filterHomologs(df, 70, 90)
        # too lenient 90 percent coverage 60-70 identity      
        
        
        
        
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
        
        print("Phrase frequencies: ")
        print(list_of_phrases_to_frequencies(temp1))
        print("________________________________") 
        print("Word frequencies: ")
        print(list_of_phrases_to_word_to_frequencies(temp1))
        print("________________________________") 
        
        
        addChemicalsToOperon(inputs_operon)
        addChemicalsToOperondf(homolog_operons)
        
        print("This is how the chemicals are added to the genes in the operon:")
        if (inputs_operon != "EMPTY"):
            for gene in inputs_operon["operon"]:
                if "chemicals" in gene.keys():
                    print(gene["chemicals"])
        print("________________________________")       
        
        print("This is how the chemicals are added to the genes in the homolog operons:")        
        if (homolog_operons is not None):
            for index, row in homolog_operons.iterrows():
                op = row['Operon']
                if op != 'EMPTY':
                    for gene in op['operon']:
                        if "chemicals" in gene.keys():
                            print(gene['chemicals'])
        print("________________________________")   
                      
        
        print("These are all the chemicals in the operons:")
        chem, chemicals_in_input_reg_operon, total_enz_num, enz_w_lig_num, total_rxn_num = getAllChemicalsInOperons(inputs_operon, homolog_operons)
        print(chem)
        print("________________________________")      
        
        print("The total number of enzymes that have chemicals found in the query operon + all homolog operons is: " + str(enz_w_lig_num))  
        print("________________________________")  
        
        print("These are all the chemicals ranked:")
        weights = {}
        weights["distance"] = 5
        weights["identity"] = 11
        weights["coverage"] = 20
        weights["frequency"] = 170
        weights["max input operon penalty"] = 250
        output_chemicals = rankChemicals(chem, chemicals_in_input_reg_operon, 0, enz_w_lig_num, weights)
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
    
    print(get_uniprot_id("NP_391277.1"))
    print("_____")
    print(get_ncbi_id("AGA23232.1"))
    print("_____")