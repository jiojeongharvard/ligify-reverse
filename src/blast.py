import subprocess
import requests
import json

from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from pprint import pprint
import pandas as pd
import streamlit as st

import xml.etree.ElementTree as ET
import re


#1) using diamond w/ database (only transcription factors) could be missing some info 2) blast of all known proteins remote-blast (5 mins)
# regulators with unknown ligands and test out in lab

    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        fasta = response.text.split("\n")
        fasta = [i for i in fasta if len(i) != 0]
        fasta = "".join(i for i in fasta if i[0] != ">")
        return fasta
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        st.error("RefSeq ID is invalid")
        return None


def uniprotID2sequence(ID: str):
    URL = f"https://rest.uniprot.org/uniprotkb/{ID}?format=json&fields=sequence"
    response = requests.get(URL)
    if response.ok:
        seq = json.loads(response.text)["sequence"]["value"]
        return seq
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        st.error("Uniprot ID is invalid")
        return None


#@st.cache_data(show_spinner=False)
def blast(acc, input_method, params, max_seqs):

    if input_method == "RefSeq":
        seq = accID2sequence(acc)
    elif input_method == "Uniprot":
        seq = uniprotID2sequence(acc)
    else:
        seq = acc
        
    flags = 'sseqid pident qcovhsp'
        # Must set this memory limit for running on a 1GB EC2 instance
    memory_limit = 0.15
  
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    log = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")
    with open(query.name, "r") as file_handle:  #opens BLAST file
        show = file_handle.readlines()
    
    # Select database to blast
    diamond_db = "/Users/jiojeong/Desktop/bHTH_RefSeq.dmnd"
    
    subprocess.call(f'diamond blastp -d {diamond_db} -q {query.name} -o {tmp.name} --outfmt 6 {flags} -b {memory_limit}'
                    f' --id {params["ident_cutoff"]} --query-cover {params["cov_cutoff"]} --max-target-seqs {max_seqs} >> {log.name} 2>&1' , shell=True)
    
    with open(tmp.name, "r") as file_handle:  #opens BLAST file
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    inDf = inDf.apply(pd.to_numeric, errors='ignore')

    try:
        inDf['sseqid'] = inDf['sseqid'].str.split("|", n=2, expand=True)[1]
    except (ValueError, KeyError):
        pass

    inDf.rename(columns= {'sseqid': 'Uniprot Id'}, inplace=True)
    inDf.rename(columns= {'pident': 'Identity'}, inplace=True)
    inDf.rename(columns= {'qcovhsp': 'Coverage'}, inplace=True)


    

    return inDf



def blast_remote(refseq_id, database='nr', num_alignments=10):
    # Perform BLAST search
    result_handle = NCBIWWW.qblast(program='blastp', database=database, sequence=refseq_id, entrez_query='txid2[Organism]', alignments=num_alignments)

    # Parse the result
    blast_records = NCBIXML.parse(result_handle)
    
    # Extract relevant information and store in a DataFrame
    results = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            print("Alignment details:")
            print("  Hit ID:", alignment.hit_id)
            print("  Hit Definition:", alignment.hit_def)
            print("  Hit Length:", alignment.length)
            print("  Number of HSPs:", len(alignment.hsps))
            for hsp in alignment.hsps:
                print("  HSP Score:", hsp.score)
                print("  HSP E-value:", hsp.expect)
                print("  HSP Query Start:", hsp.query_start)
                print("  HSP Query End:", hsp.query_end)
                print("  HSP Hit Start:", hsp.sbjct_start)
                print("  HSP Hit End:", hsp.sbjct_end)
                print("  HSP Alignment Length:", hsp.align_length)
                print("  HSP Identity:", hsp.identities)
                print("  HSP Query Sequence:", hsp.query)
                print("  HSP Hit Sequence:", hsp.sbjct)
                print()
                
                # Calculate identity and coverage percentages
                
                alignment_length = hsp.align_length
                identity = hsp.identities
                coverage = (identity / alignment_length) * 100

                # Extract UniProt ID from subject ID if it follows the UniProt format
                subject_id = None
                if '|' in alignment.hit_id:
                    if alignment.hit_id.split('|')[0] == "ref": 
                        subject_id = alignment.hit_id.split('|')[1]

                results.append({
                    'NCBI Id': subject_id,
                    'Identity': (identity / alignment_length) * 100,
                    'Coverage': coverage
                })

    df = pd.DataFrame(results)
    df_filtered = df[df['NCBI Id'].notnull()]
    return df_filtered


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
 
if __name__ == "__main__":
        
    # # Example usage
    # refseq_id = "WP_004925500.1"
    # blast_df = blast_remote(refseq_id)
    # print(blast_df)
        
    print(get_ncbi_id("WP_004925500.1"))

    # if acc != None:
    #     df = blast(acc)
    #     pprint(df)
    # else:
    #     print("bad seq")

    # print(seq)