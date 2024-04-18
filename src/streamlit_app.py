import streamlit as st

import re
import pandas as pd
import requests
import json
from chemical import get_smiles_from_chebi, smiles_to_image
from operon import getOperon, filterRegFromOperon, filterRegFromOperondf, get_enzyme_description, get_enzyme_description_df, addChemicalsToOperon, addChemicalsToOperondf, getAllChemicalsInOperons, rankChemicals, list_of_phrases_to_frequencies
from blast import blast, blast_remote
from accID2operon import acc2operon
import datetime
from concurrent.futures import ThreadPoolExecutor, TimeoutError

st.set_page_config(page_title="Ligify", layout='wide', initial_sidebar_state='auto', 
    #page_icon="images/Snowprint_favicon.png"
    )

hide_streamlit_style = '''
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
'''
st.markdown(hide_streamlit_style, unsafe_allow_html=True)


# Removes border around forms
css = r'''
    <style>
        [data-testid="stForm"] {border: 0px}
    </style>
'''
st.markdown(css, unsafe_allow_html=True)


# Removes the full-screen button for various elements
style_fullscreen_button_css = """
    button[title="View fullscreen"] {
        display: none;
    }
    button[title="View fullscreen"]:hover {
        display: none;
        }
    """
st.markdown(
    "<style>"
    + style_fullscreen_button_css
    + "</styles>",
    unsafe_allow_html=True,
)



# Initialize state variables
if "data" not in st.session_state:
        st.session_state.data = False

if 'SUBMITTED' not in st.session_state:
    st.session_state.SUBMITTED =  False


def _connect_form_cb(connect_status):
    st.session_state.SUBMITTED = connect_status
    st.session_state.data = False



# HEADER: Title and basic input
head = st.container()
head1, head2, head3 = head.columns((1,2,1))

# head2.image("images/Snowprint_Logo.png", use_column_width=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Predict a regulator's inducer molecule</h3>", unsafe_allow_html=True)

selection_container = st.container()
sel1, sel2, sel3 = selection_container.columns((1,2,1))

input_method = sel2.radio(label="Choose an input format", \
        options=("RefSeq", "Uniprot", "Protein sequence"))

input_container = st.container()
in1, in2, in3 = input_container.columns((1,2,1))

if input_method == "RefSeq":
    acc = in2.text_input(label="RefSeq ID", value="WP_013083972.1", label_visibility="hidden")
elif input_method == "Uniprot":
    acc = in2.text_input(label="UniprotID", value="P43506", label_visibility="hidden")
elif input_method == "Protein sequence":
    protein_seq_input = in2.text_area(label="Protein sequence", height=200, label_visibility="hidden")
    if re.match(r'^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]*$', protein_seq_input) and len(protein_seq_input) > 50:
        acc = protein_seq_input
    else:
        in2.error("Protein sequence must only contain amino acid characters and have a length over 50. Make sure there is no new line after the sequence.")




# # Advanced options
# options = st.container()
# options1, options2, options3 = options.columns((1,3,1))

# #with options2:
    


# Side bar
with st.sidebar:
    st.write("Prokaryotic transcription factors can be repurposed as chemical measurement tools for synthetic biology.")
    
    st.write("To repurpose a transcription factor, the specific inducer molecule it binds to must be determined.")

    st.write("Ligify predicts transcription factor-ligand interactions by analyzing conserved enzymes in local genomic contexts.")

    # GitHub and Email links
    st.markdown("<p style='font-size: 12px'>If you have any questions or would like to report any bugs, please contact us via <a href='mailto: simonsnitz@gmail.com'>Email</a>. \
        Our code is publically available on <a href='https://github.com/simonsnitz/snowstream'>GitHub</a>.</p>", unsafe_allow_html=True)

    # st.markdown("<div style='font-size: 12px;'>d'Oelsnitz S., Stofel S.K., and Ellington A.D. (2023) Snowprint: a predictive tool for genetic biosensor discovery. \
    #             <i>bioRxiv</i> <b>DOI:</b><a href='https://www.biorxiv.org/content/10.1101/2023.04.29.538814v1'>10.1101/2023.04.29.538814v1</a></div> <br>", unsafe_allow_html=True)

    st.divider()
    
    # Display a header or some text indicating these are advanced options
    st.header("Advanced Options")
    st.write("Run Mode")
    
    blast_mode = st.selectbox("Would you like to search for homologs? If so, would you like to run using a DIAMOND database (faster) or remote BLAST (slower but more homologs)?",
   ("NO homolog search",  "Homolog search: DIAMOND database", "Homolog search: remote BLAST"), 
   index=0,
   placeholder="Select method...",
)
    st.divider()
    
    st.write("BLAST Parameters (Ignore if not running homology search)")
    
    
        
    # Blast Parameters with default values
    ident_cutoff = st.slider("Identity Cutoff (%)", min_value=0, max_value=100, value=70, key='ident_cutoff')
    cov_cutoff = st.slider("Coverage Cutoff (%)", min_value=0, max_value=100, value=90, key='cov_cutoff')

    # Checkbox for filtering redundant results with a default value
    filter_redundant = st.checkbox("Filter Redundant Homologs", value=True, key='filter_redundant')

    # Number input for maximum number of homologs with a default value
    max_homologs = st.number_input("Maximum Number of Homologs", min_value=1, max_value=100, value=50, key='max_homologs')
    st.divider()
    
    st.write("Chemical Scoring")
    
    st.markdown("<p style='font-size: 10px'>The ligands excluded from ranking are the following: 1, 5-dihydroflavin, acetyl-CoA(4-), acyl-CoA(4-), ADP(3-), aldehyde, ammonium, AMP(2-), AMP 3-end(1-) residue, ATP(4-), carbon dioxide, carboxylic acid anion, coenzyme A(4-), copper(1+), copper(2+), dioxygen, diphosphate(3-), DNA 5-phosphate polyanion, FAD(3-), FADH(2)(2-), flavin(1-), FMN(3-), FMNH(2)(2-), GDP(3-), GMP(2-), GTP(4-), hydrogen acceptor, hydrogen donor, hydrogen peroxide, hydrogenphosphate, H group, hydron, iron(2+), iron(3+), lipid II(3-), NAD(1-), NADH(2-), NADP(3-), NADPH(4-), potassium(1+), sodium(1+), superoxide, thiol group, water</p>", unsafe_allow_html=True)
      
    
    chemical_score_cutoff = st.slider("Score Cutoff", min_value=0, max_value=100, value=50, key='chemical_score_cutoff')
    
    distance_weight = st.slider("Weight: Distance From Regulator in Operon", min_value=0, max_value=50, value=5, key='distance_weight')
    identity_weight = st.slider("Weight: Identity (%)", min_value=0, max_value=50, value=11, key='identity_weight')
    coverage_weight = st.slider("Weight: Coverage (%)", min_value=0, max_value=50, value=20, key='coverage_weight')
    num_occurrence_weight = st.slider("Weight: Number of Appearances", min_value=0, max_value=50, value=17, key='num_occurrence_weight')
    num_occurrence_weight = num_occurrence_weight * 10
    query_operon_max_penalty = st.slider("Penalty for Missing in Query Protein's Operon", min_value=0, max_value=50, value=25, key='query_operon_max_penalty')
    query_operon_max_penalty = query_operon_max_penalty * 10
    
    

    

weights = {"distance": distance_weight, "identity": identity_weight, "coverage": coverage_weight, "frequency": num_occurrence_weight, "max input operon penalty": query_operon_max_penalty}
    
blast_params = {"ident_cutoff": ident_cutoff, "cov_cutoff": cov_cutoff}

search_homolog = False
use_blast_remote = False

if (blast_mode == "NO homolog search"):
    search_homolog = False
    use_blast_remote = False
elif (blast_mode == "Homolog search: remote BLAST"):
    search_homolog = True
    use_blast_remote = True
elif (blast_mode == "Homolog search: DIAMOND database"):
    search_homolog = True
    use_blast_remote = False
    

# FORM
with st.form(key='ligify'):

    # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([4,1,4])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))


# RUN LIGIFY
if st.session_state.SUBMITTED:
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    st.write(f"Submitted at {current_time}")
    
    chemical = st.container()
    chem1, chem2, chem3 = chemical.columns((2,2,2))
    chemical_col, spacer = chemical.columns((4,1))
    
    intermediate = st.container()
    input1, input2 = intermediate.columns((4,1))
    blast_col, enzyme_col = intermediate.columns((3,3))
    
    metrics = st.container()
    m_spacer1, metrics_col, m_spacer2 = metrics.columns((1,7,1))
    
    if (search_homolog and input_method == "Protein sequence"):
        st.error("When input format is protein sequence, must perform homology search")

    with st.spinner("Running..."):
        
        prog_container = st.container()
        prog_spacerL, prog, prog_spacerR = prog_container.columns((1,3,1))
        prog_bar = prog.progress(10, text="1. Fetching homologs")
        blast_df = pd.DataFrame()
        
        def execute_blast():
            if search_homolog and not use_blast_remote:
                return blast(acc, input_method, blast_params, 100)
            elif search_homolog and use_blast_remote:
                return blast_remote(acc, input_method, blast_params, 100)
            else:
                return pd.DataFrame()  # or appropriate default action

        # Using ThreadPoolExecutor to run blast functions with a timeout
        with ThreadPoolExecutor() as executor:
            future = executor.submit(execute_blast)
            try:
                # Set timeout to 3600 seconds (1 hour)
                blast_df = future.result(timeout=3600)
            except TimeoutError:
                st.error("The BLAST operation timed out after an hour. Please try again later.")
            
        # If BLAST does not return anything, troubleshoot the issue.
        if blast_df.empty and search_homolog:
            if input_method == "Protein sequence":
                blast_col.error("BLAST failed. Try running a RefSeq or Uniprot ID for more detailed error codes")
            else:
                # mode, genes = troubleshoot(input_method, acc)
                mode = "genes"
                if mode == "not bacteria":
                    st.error("Protein is not from Bacteria. Ligify only works for bacterial proteins")
                elif mode == "genes":
                    st.error("No blast results returned. Try running the command line tool https://github.com/simonsnitz/Snowprint")

        # BLAST results look good aka not empty
        else:

            #if 'filter redundant' box checked, filter out homologs that have the same %identity and %coverage
            def filter_blastDf(blast_df):
                homolog_dict = []
                ident_covs = []
                for i, row in blast_df.iterrows():
                    if (use_blast_remote):
                        entry =   {"NCBI Id": row["NCBI Id"], "Identity": row["Identity"],"Coverage": row["Coverage"]}
                    else:
                        entry =   {"Uniprot Id": row["Uniprot Id"], "Identity": row["Identity"],"Coverage": row["Coverage"]}
                    to_compare =   {"Identity": row["Identity"],"Coverage": row["Coverage"]}
                    if to_compare not in ident_covs:
                        homolog_dict.append(entry)
                        ident_covs.append(to_compare)
                return homolog_dict
            
            prog_bar.progress(20, text=f"2. Filtering redundant homologs")
            
 
            if (search_homolog):
                if filter_redundant:
                    homolog_dict = filter_blastDf(blast_df)
                else:
                    if (use_blast_remote):
                        homolog_dict = [
                            {"NCBI Id": row["NCBI Id"], "Identity": row["Identity"],"Coverage": row["Coverage"]}
                            for i, row in blast_df.iterrows()
                        ]
                    else:
                        homolog_dict = [
                            {"Uniprot Id": row["Uniprot Id"], "Identity": row["Identity"],"Coverage": row["Coverage"]}
                            for i, row in blast_df.iterrows()
                        ]
            # limit search to specified number of homologs
            if (search_homolog):
                if (len(homolog_dict) > max_homologs):
                    homolog_dict = homolog_dict[0:max_homologs]
            
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


            # Create an info section on the input protein
            def uniprotID2info(ID: str):
                URL = f"https://rest.uniprot.org/uniprotkb/{ID}?format=json&fields=sequence,organism_name,protein_name"
                response = requests.get(URL)
                if response.ok:
                    data = json.loads(response.text)
                    out = {
                        "Annotation": data["proteinDescription"]["recommendedName"]["fullName"]["value"],
                        "Organism": data["organism"]["scientificName"],
                        "Lineage": data["organism"]["lineage"],
                    }
                    return out
                else:
                    print("FATAL: Bad uniprot API request "+ str(response.status_code))
                    st.error("Uniprot ID is invalid")
                    return None
                
            prog_bar.progress(30, text=f"3. Fetching information on input protein")
            uniprot = ""
            if (input_method == "RefSeq"):
                uniprot = get_uniprot_id(acc)
            elif (input_method == "Uniprot"):
                uniprot = acc
            else:
                if (search_homolog):
                    if (blast_df.iloc[0]['Coverage'] == 100 and blast_df.iloc[0]['Identity'] == 100):
                        if (not use_blast_remote):
                            uniprot = blast_df.iloc[0]["Uniprot Id"]
                        else:
                            uniprot = get_uniprot_id(blast_df.iloc[0]["NCBI Id"])
            
            if (input_method == "Protein sequence"):
                acc = uniprot
            
            try:
                protein_data = uniprotID2info(uniprot)
                input1.subheader("Input")
                input1.markdown("Annotation: "+protein_data["Annotation"])
                input1.markdown("Organism: "+protein_data["Organism"])
                lineage = "".join(i+", " for i in protein_data["Lineage"])[:-2]
                input1.markdown("Lineage: "+lineage)
            except:
                pass

            blast_col.subheader("BLAST results")
            
            num_homologs = 0
            if (search_homolog):
                homolog_df = pd.DataFrame(homolog_dict)
                homolog_df = homolog_df[~((homolog_df['Identity'] == 100) & (homolog_df['Coverage'] == 100))]
                num_homologs = homolog_df.shape[0]

                blast_col.dataframe(homolog_df)
            
            inputs_operon = "EMPTY"
            
            prog_bar.progress(40, text=f"4. Fetching operons")
            if (search_homolog):
                inputs_operon, homolog_operons = getOperon(acc, homolog_df)
            else:
                inputs_operon = acc2operon(acc)
                homolog_operons = None
               
            if (inputs_operon == "EMPTY"):
                st.error("Operon of input protein cannot be found.")
            
            prog_bar.progress(50, text=f"5. Filtering regulators out of operons")
            filterRegFromOperon(acc, inputs_operon)
            
            if (search_homolog):
                filterRegFromOperondf(homolog_operons)   
            
            prog_bar.progress(60, text=f"6. Fetching enzyme descriptions")    
            enzyme_descriptions = get_enzyme_description(inputs_operon)
            
            prog_bar.progress(70, text=f"7. Adding chemicals to operons")    
            addChemicalsToOperon(inputs_operon)
            
            if (search_homolog):
                temp = get_enzyme_description_df(homolog_operons)
                enzyme_descriptions.extend(temp)
                addChemicalsToOperondf(homolog_operons)
            
            enzyme_description_frequencies = list_of_phrases_to_frequencies(enzyme_descriptions)
            
            chem, chemicals_in_input_reg_operon, total_enz_num, enz_w_lig_num, total_rxn_num  = getAllChemicalsInOperons(inputs_operon, homolog_operons)
            
            prog_bar.progress(80, text=f"8. Ranking chemicals")  
            output_chemicals = rankChemicals(chem, chemicals_in_input_reg_operon, 0, enz_w_lig_num, weights)
            num_chemicals = len(output_chemicals)
            
            # need name, structure, score
            
            if (len(output_chemicals) >= 1):
                chem1.subheader("1st place")
                chem1_img = smiles_to_image(get_smiles_from_chebi(output_chemicals[0]["ChEBI ID"]))
                chem1.image(chem1_img, caption=f"name: {output_chemicals[0]['Chemical Name']}, score: {output_chemicals[0]['Score']}")
                
            if (len(output_chemicals) >= 2):
                chem2.subheader("2nd place")
                chem2_img = smiles_to_image(get_smiles_from_chebi(output_chemicals[1]["ChEBI ID"]))
                chem2.image(chem2_img, caption=f"name: {output_chemicals[1]['Chemical Name']}, score: {output_chemicals[1]['Score']}")
            
            if (len(output_chemicals) >= 3):
                chem3.subheader("3rd place")
                chem3_img = smiles_to_image(get_smiles_from_chebi(output_chemicals[2]["ChEBI ID"]))
                chem3.image(chem3_img, caption=f"name: {output_chemicals[2]['Chemical Name']}, score: {output_chemicals[2]['Score']}")
            
            
            chemical_col.subheader("Chemical scoring")
            output_chemicals_df = pd.DataFrame(output_chemicals)
            if ('Subscore' in output_chemicals_df.columns):
                output_chemicals_df = output_chemicals_df.drop(columns=['Subscore'])
            chemical_col.dataframe(output_chemicals_df)
            
            # # Creating the bar chart
            # fig, ax = plt.subplots()
            # ax.bar(enzyme_description_frequencies.keys(), enzyme_description_frequencies.values())
            # ax.set_xlabel('Words')
            # ax.set_ylabel('Frequency')
            # ax.set_title('Word Frequency Diagram')
            
            enzyme_col.subheader("Summary of all the enzymes found")
            # Display the plot in the right column of the Streamlit app
            
            enzyme_df = pd.DataFrame(list(enzyme_description_frequencies.items()), columns=['Enzyme', 'Frequency'])
    
    
            enzyme_col.dataframe(enzyme_df)
            
            prog_bar.empty()
            
            metrics_col.subheader("Search metrics")
            m_homologs, m_tot_genes, m_reactions, m_chemicals = metrics_col.columns(4)
            m_homologs.metric("Homologs", num_homologs)
            m_tot_genes.metric("Total proteins", total_enz_num)
            m_reactions.metric("Reactions found", total_rxn_num)
            m_chemicals.metric("Unique chemicals", num_chemicals)
            metrics_col.divider()
                
            
            
            
            
            