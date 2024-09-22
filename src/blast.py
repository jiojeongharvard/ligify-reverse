import subprocess
import requests
import json

from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd
import streamlit as st

from Bio.Blast.Applications import NcbiblastpCommandline

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
    
    if (seq is None):
        return pd.DataFrame()
        
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
    diamond_db = "/Users/jiojeong/Documents/ligify-reverse-1/bHTH_RefSeq.dmnd"
    
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


def blast_remote(acc: str, input_method, params, num_aligns):

    # print("NOTE: Starting BLAST")
    
    if input_method == "RefSeq":
        seq = accID2sequence(acc)
    elif input_method == "Uniprot":
        seq = uniprotID2sequence(acc)
    else:
        seq = acc
    
    if (seq is None):
        return pd.DataFrame()

    blast_db = "nr"

    if seq != None:
            # Must have BLAST+ executables in PATH to run this
        blast_cline = NcbiblastpCommandline(db=blast_db, outfmt="6 sseqid pident qcovs", \
            num_alignments=num_aligns, remote=True)

        results, err = blast_cline(stdin=seq)

        results = results.split("\n")[:-1]
        
        homologs = [{"NCBI Id": r.split("|")[1], \
                "Identity": r.split("|")[2].split("\t")[1], \
                "Coverage": r.split("|")[2].split("\t")[2].strip()} \
                for r in results ]

        df = pd.DataFrame(homologs)
        df_filtered = df[df['NCBI Id'].notnull()]
        df_filtered_end = df_filtered[df_filtered['NCBI Id'].str.len() != 4]
        df_filtered_end['Identity'] = df_filtered_end['Identity'].astype(float)
        df_filtered_end['Coverage'] = df_filtered_end['Coverage'].astype(float)
        final_df = df_filtered_end[(df_filtered_end['Identity'] > params["ident_cutoff"]) & (df_filtered_end['Coverage'] > params["cov_cutoff"])]
        return final_df

    else:
        return None
    

# def blast_remote_old(refseq_id, num_alignments, database='nr'):
#     # Perform BLAST search
#     seq = accID2sequence(refseq_id)
#     result_handle = NCBIWWW.qblast(program='blastp', database=database, sequence=seq, entrez_query='txid2[Organism]', alignments=num_alignments)

#     # Parse the result
#     blast_records = NCBIXML.parse(result_handle)
    
#     # Extract relevant information and store in a DataFrame
#     results = []
#     for blast_record in blast_records:
#         for alignment in blast_record.alignments:
#             print("Alignment details:")
#             print("  Hit ID:", alignment.hit_id)
#             print("  Hit Definition:", alignment.hit_def)
#             print("  Hit Length:", alignment.length)
#             print("  Number of HSPs:", len(alignment.hsps))
#             for hsp in alignment.hsps:
#                 print("  HSP Score:", hsp.score)
#                 print("  HSP E-value:", hsp.expect)
#                 print("  HSP Query Start:", hsp.query_start)
#                 print("  HSP Query End:", hsp.query_end)
#                 print("  HSP Hit Start:", hsp.sbjct_start)
#                 print("  HSP Hit End:", hsp.sbjct_end)
#                 print("  HSP Alignment Length:", hsp.align_length)
#                 print("  HSP Identity:", hsp.identities)
#                 print("  HSP Query Sequence:", hsp.query)
#                 print("  HSP Hit Sequence:", hsp.sbjct)
#                 print()
                
#                 # Calculate identity and coverage percentages
                
                

#                 # Extract UniProt ID from subject ID if it follows the UniProt format
#                 subject_id = None
#                 if '|' in alignment.hit_id:
#                     if alignment.hit_id.split('|')[0] == "ref": 
#                         subject_id = alignment.hit_id.split('|')[1]

#                 results.append({
#                     'NCBI Id': subject_id,
#                     'Identity': (hsp.identities / hsp.align_length) * 100,
#                     'Coverage': (hsp.align_length / len(seq)) * 100
#                 })

#     df = pd.DataFrame(results)
#     df_filtered = df[df['NCBI Id'].notnull()]
#     # take out the pdb ones
#     return df_filtered


if __name__ == "__main__":
        
    # Example usage
    refseq_id = "WP_004925500.1"
    print(accID2sequence(refseq_id))
        
