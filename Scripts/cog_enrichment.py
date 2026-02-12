import pandas as pd
from tqdm import tqdm
from scipy.stats import fisher_exact


"""
This script processess the results from eggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) and
performed an enrichement analysis of COG categories between genomes with or without T6SS for a list of family hosts.

Input parameters:
    -Path to the files containing eggNOG-mapper results for chromosomes/plasmids with or without T6SS:
    {family}_{genome}_{t6ss}_eggnogmapper.tsv

Outputs:
    -.csv file with frequency of COG categories for T6SS+ and T6SS- genomes (one for each family hots and genome type)
    -.csv file with enrichment results (one for each family hots and genome type)
"""


def assign_t6ss_comp(row):
    """This function assign proteins as T6SS components to remove them from the enrichment analysis"""
    description = row["Description"]
    if "type VI" in description or "Type VI" in description:
        if description =="LVIVD repeat":
            valor = 0
        else:
            valor = 1
    elif "FHA" in description:
        valor = 1
    else:
        valor = 0
    return valor

def process_eggnogmapper(filename):
    """This function process annotation results using Eggnogmapper"""
   
    dicc_catego={"J": "Translation, Ribosomal Structure and Biogenesis", "A": "RNA Processing and Modification", "K": "Transcription", "L": "Replication, Recombination and Repair", "B":"Chromatin Structure and Dynamics", "D": "Cell Cycle Control, Cell Division, Chromosome Partitioning", "Y": "Nuclear structure", "V": "Defense Mechanisms", "T": "Signal Transduction Mechanisms", "M": "Cell Wall/Membrane/Envelope Biogenesis", "N": "Cell Motility", "Z": "Cytoskeleton", "W": "Extracellular Structures", "U": "Intracellular trafficking, Secretion and Vesicular Transport", "O":"Posttranslational Modification, Protein Turnover, Chaperones", "X":"Mobilome: Prophages, Transposons", "C":"Energy Production and Conversion", "G": "Carbohydrate Transport and Metabolism", "E": "Amino Acid Transport and Metabolim", "F": "Nucleotide Transport and Metabolism", "H":"Coenzyme Transport and Metabolism", "I": "Lipid Transport and Metabolism", "P": "Inorganic Ion Transport and Metabolism", "Q": "Secondary Metabolites Biosynthesis, Transport and Catabolism", "R": "General Function Prediction Only", "S": "Function Unknown", "-": "-"}
    
    df1 = pd.read_csv(filename, sep="\t")
    df1["COG_catego_desc"]=df1["COG_category"].apply(lambda x: dicc_catego[x] if x in dicc_catego else "-")
    tqdm.pandas(desc="my bar!")
    df1["T6SS_Component"]=df1.progress_apply(assign_t6ss_comp, axis=1)
    # remove T6SS components from analysis
    df1_filt = df1[df1["T6SS_Component"]!=1]
    list_catego = df1_filt["COG_category"].to_list()
    # if there are proteins associated with >1 biological process/functions, count them separate
    dicc_aux={}
    for catego in list_catego:
        if catego in dicc_catego:
            if catego in dicc_aux:
                dicc_aux[catego]+=1
            else:
                dicc_aux[catego]=1
        else:
            for i in catego:
                if i in dicc_catego:
                    if i in dicc_aux:
                        dicc_aux[i]+=1
                    else:
                        dicc_aux[i]=1           

    # build dataframe with individual COG categories
    lista_dicc = []
    for catego in dicc_aux:
        lista_dicc.append({
            "COG_category" : catego,
            "COG_category_desc" : dicc_catego[catego],
            "Count" : dicc_aux[catego]
        })
    
    df_freq = pd.DataFrame(lista_dicc)
    return df_freq


def contingency_table(df_t6ss, df_not6ss, var_family, genome):
    """This function performs an enrichment analysis of COG categories between two groups"""
    
    OUTPUT1 = f"{var_family}_{genome}_T6SS_vs_noT6SS_COG_categories_TABLE.tsv"
    OUTPUT2 = f"{var_family}_{genome}_T6SS_vs_noT6SS_COG_categories_enrichment_ODDSRATIO.tsv"

    dicc_catego={"J": "Translation, Ribosomal Structure and Biogenesis", "A": "RNA Processing and Modification", "K": "Transcription", "L": "Replication, Recombination and Repair", "B":"Chromatin Structure and Dynamics", "D": "Cell Cycle Control, Cell Division, Chromosome Partitioning", "Y": "Nuclear structure", "V": "Defense Mechanisms", "T": "Signal Transduction Mechanisms", "M": "Cell Wall/Membrane/Envelope Biogenesis", "N": "Cell Motility", "Z": "Cytoskeleton", "W": "Extracellular Structures", "U": "Intracellular trafficking, Secretion and Vesicular Transport", "O":"Posttranslational Modification, Protein Turnover, Chaperones", "X":"Mobilome: Prophages, Transposons", "C":"Energy Production and Conversion", "G": "Carbohydrate Transport and Metabolism", "E": "Amino Acid Transport and Metabolim", "F": "Nucleotide Transport and Metabolism", "H":"Coenzyme Transport and Metabolism", "I": "Lipid Transport and Metabolism", "P": "Inorganic Ion Transport and Metabolism", "Q": "Secondary Metabolites Biosynthesis, Transport and Catabolism", "R": "General Function Prediction Only", "S": "Function Unknown", "-": "-"}

    # build dataframe T6SS+ vs T6SS-
    df = pd.DataFrame(0, index=dicc_catego.keys(), columns=["T6SS-", "T6SS+"])
    categos_t6ss = df_t6ss["COG_category"].to_list()
    categos_not6ss = df_not6ss["COG_category"].to_list()
    
    for cog in dicc_catego:
        count_not6ss = 0
        count_t6ss = 0
        if cog in categos_t6ss:
            count_t6ss = int(df_t6ss[df_t6ss["COG_category"]==cog]["Count"].to_list()[0])
            
        if cog in categos_not6ss:
            count_not6ss = int(df_not6ss[df_not6ss["COG_category"]==cog]["Count"].to_list()[0])
           
        df.at[cog, "T6SS-"]=count_not6ss
        df.at[cog, "T6SS+"]=count_t6ss

    df.to_csv(OUTPUT1, sep="\t")
        
    # filter null values in both columns
    df = df[(df["T6SS+"] > 0) & (df["T6SS-"] > 0)]

    ##########################################
    ##########################################
    ##### FUNCTIONAL ENRICHMENT ANALYSIS #####
    ##########################################
    ##########################################

    results = []
    for index, row in df.iterrows():
        # build 2x2 contingency table for each category
        table = [[row["T6SS+"], df["T6SS+"].sum() - row["T6SS+"]],
                [row["T6SS-"], df["T6SS-"].sum() - row["T6SS-"]]]
              
        # Fisher exact test
        # An odds ratio greater than 1 indicates that the condition or event is more likely to occur in the first group
        odds_ratio, p_valor_fisher = fisher_exact(table, alternative="greater")
        
        results.append({
            "COG_category": index,
            "Odds_ratio": odds_ratio,
            "p-value Fisher": p_valor_fisher
        })


    results_df = pd.DataFrame(results)
    # Benjamini-Hochberg correction
    results_df["Rank Fisher"] = results_df["p-value Fisher"].rank(method="min")
    results_df["p-adj Fisher BH"] = results_df["p-value Fisher"] * len(results_df) / results_df["Rank Fisher"]
    results_df["p-adj Fisher BH"] = results_df["p-adj Fisher BH"].clip(upper=1)  
    results_df["p-adj Fisher BH significant"] = results_df["p-adj Fisher BH"].apply(lambda x: "Yes" if x <0.05 else "No")
    results_df.to_csv(OUTPUT2, sep="\t", index=False)
    results_signif = results_df[results_df["p-adj Fisher BH"] < 0.05]

    return results_signif



list_families = ["Campylobacteraceae", "Yersiniaceae", "Vibrionaceae", "Rhodobacteraceae", "Burkholderiaceae", "Enterobacteriaceae"]
for family in list_families:
    # process eggNOG-mapper raw results
    df_chr_not6ss = process_eggnogmapper(f"{family}_CHR_noT6SS_eggnogmapper.tsv", family, "CHR_noT6SS")
    df_chr_t6ss = process_eggnogmapper(f"{family}_CHR_T6SS_eggnogmapper.tsv", family, "CHR_T6SS")
    df_pla_not6ss = process_eggnogmapper(f"{family}_PLA_noT6SS_eggnogmapper.tsv", family, "PLA_noT6SS")
    df_pla_t6ss = process_eggnogmapper(f"{family}_PLA_T6SS_eggnogmapper.tsv", family, "PLA_T6SS")

    # enrichment analysis for each genomic platform (chromosome, plasmid)
    df_enrich_chr = contingency_table(df_chr_t6ss, df_chr_not6ss, family, "CHR")
    df_enrich_pla = contingency_table(df_pla_t6ss, df_pla_not6ss, family, "PLA")

    # comparison of enriched COG for each genomic platform 
    cogs_chr = df_enrich_chr["COG_category"].to_list()
    cogs_pla = df_enrich_pla["COG_category"].to_list()
    only_chr = set(cogs_chr) - set(cogs_pla)
    only_pla = set(cogs_pla) - set(cogs_chr)
    common = set(cogs_chr) & set(cogs_pla)
    print(f"\nCOG categories enriched in T6SS+ genomes for {family}: ")
    print("CHR: ", len(only_chr))
    print(only_chr)
    print("\nPLASMIDS: ", len(only_pla))
    print(only_pla)
    print("\nCHR and PLASMIDS: ", len(common))
    print(common)

print("\nFINISHED!")


