import streamlit as st
from src.blast import blast
# from src.get_genome_coordinates import get_genome_coordinates, get_genome_coordinates_batch
# from src.accID2operon import acc2operon
# from src.fetch_promoter import fetch_promoter
# from src.fetch_operator import fetch_operator
# from src.troubleshoot import troubleshoot

import re
import pandas as pd
import requests
import json


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
        in2.error("Protein sequence must only contain amino acid characters and have a length over 50")




# Advanced options
options = st.container()
options1, options2, options3 = options.columns((1,3,1))



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



# Format advanced options
blast_params = {
    "ident_cutoff": 60,
    "cov_cutoff": 95
}
filter_redundant = True
max_homologs = 30



# FORM
with st.form(key='ligify'):

        # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))


# RUN LIGIFY
if st.session_state.SUBMITTED:
    st.write("submitted")

    intermediate = st.container()
    input1, input2 = intermediate.columns((4,1))
    blast_col, coordinates_col, homologs_col = intermediate.columns((1.3,2,2.2))

    with st.spinner("blasting your protein"):


        blast_df = blast(acc, input_method, blast_params, max_seqs=500)

        # If BLAST does not return anything, troubleshoot the issue.
        if blast_df.empty:
            if input_method == "Protein sequence":
                blast_col.error("BLAST failed. Try running a RefSeq or Uniprot ID for more detailed error codes")
            else:
                mode, genes = troubleshoot(input_method, acc)
                if mode == "not bacteria":
                    st.error("Protein is not from Bacteria. Snowprint only works for bacterial proteins")
                elif mode == "genes":
                    st.error("No blast results returned. Try running the command line tool https://github.com/simonsnitz/Snowprint")

        # BLAST results look good
        else:

            #if 'filter redundant' box checked, filter out homologs that have the same %identity and %coverage
            def filter_blastDf(blast_df):
                homolog_dict = []
                ident_covs = []
                for i, row in blast_df.iterrows():
                    entry =   {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
                    to_compare =   {"identity": row["Identity"],"coverage": row["Coverage"]}
                    if to_compare not in ident_covs:
                        homolog_dict.append(entry)
                        ident_covs.append(to_compare)
                return homolog_dict

            if filter_redundant:
                homolog_dict = filter_blastDf(blast_df)
            else:
                homolog_dict = [
                    {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
                    for i, row in blast_df.iterrows()
                    ]

            # limit search to specified number of homologs
            homolog_dict = homolog_dict[0:max_homologs]



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

            uniprot = blast_df.iloc[0]["Uniprot Id"]
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
            blast_col.dataframe(homolog_dict)