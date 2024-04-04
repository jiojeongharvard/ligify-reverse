# ligify-reverse

Ligify-reverse is a program that predicts protein-ligand interactions. 
It is based off on the program Ligify, which can be read more about [here]().
Ligify-reverse takes a protein (refseq id, uniprot id, or primary sequence) as an input and outputs the likely chemical that it binds to.
The user can run the program with the following command: streamlit run streamlit_app.py
The user has the option to run the program on three modes: no homology search of input protein, blast search through local diamond database, and remote blast search.

If you want to run homology search using a diamond database, you can download the diamond database [here](https://www.dropbox.com/scl/fi/4xqaymf7oxz4cqj6knage/bHTH_RefSeq.dmnd?rlkey=dhcujdxgfqlkmi6rd6oy8fngj&dl=1).
After downloading the database, please set the "diamond_db" path in blast.py to where the file is saved in your device.

The user also has the option to set the BLAST parameters for homology searching and the weights for the chemical scoring function. 

The output shown on the page will consist of the list of homologs (if searched), list of chemicals ranked by a scoring scheme, and list of all enzyme descriptions found in the operons of query protein + homologs. 





