from operon import filterHomologs, getOperon, filterRegFromOperon, filterRegFromOperondf, get_enzyme_description, get_enzyme_description_df, addChemicalsToOperon, addChemicalsToOperondf, getAllChemicalsInOperons, rankChemicals
from blast import blast
from accID2operon import acc2operon

if __name__ == "__main__":
    query = ""
    ident_cutoff = 70
    cov_cutoff = 90
    max_num_homologs = 20
    chemical_score_cutoff = 30
    
    distance_weight = 5
    identity_weight = 11
    coverage_weight = 20
    num_occurrence_weight = 170
    query_operon_max_penalty = 250
    
    weights = {}
    weights["distance"] = distance_weight
    weights["identity"] = identity_weight
    weights["coverage"] = coverage_weight
    weights["frequency"] = num_occurrence_weight
    weights["max input operon penalty"] = query_operon_max_penalty
    
    params = {"ident_cutoff": ident_cutoff, "cov_cutoff": cov_cutoff}
    
    search_homolog = True
    blast_remove = False
    
    
    inputs_operon = None
    homolog_operons = None
    if (search_homolog):
        df = blast(query, "RefSeq", params, max_num_homologs)
        filtered_df = filterHomologs(df, 70, 90)  
        inputs_operon, homolog_operons = getOperon(query, filtered_df)
    else:
        inputs_operon = acc2operon(query)
    
    filterRegFromOperon(query, inputs_operon)
    if (search_homolog):
        filterRegFromOperondf(homolog_operons)
    
    enzyme_descriptions = get_enzyme_description(inputs_operon)
    addChemicalsToOperon(inputs_operon)
    if (search_homolog):
       temp = get_enzyme_description_df(homolog_operons)
       enzyme_descriptions.extend(temp)
       addChemicalsToOperondf(homolog_operons)
    
    chem, chemicals_in_input_reg_operon, num_genes = getAllChemicalsInOperons(inputs_operon, homolog_operons)
    output_chemicals = rankChemicals(chem, chemicals_in_input_reg_operon, 0, num_genes, weights)
    print(output_chemicals)