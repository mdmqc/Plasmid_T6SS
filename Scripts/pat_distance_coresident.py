import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt


"""
This script extracts the patristic distances in the TssC tree for co-resident chromosomal and plasmid T6SSs
and represents the Cumulative Distribution Function (CDF) of those distances in comparison with the closest plasmid homolog. 
Input parameters:
    -file1: Path to the file containing TssC clusterization 99id, 100cov representatives metadata
    -file2: Path to the file containing the patristic distance matrix
    -file3: Path to the file containing the co-resident chromosomes and plasmids metadata
    -file4: Path to the file containing the closest plasmid TssC homologs to chromosomal TssC (generated with patristic_distances.py)
Outputs:
    -OUTPUT1: Co-resident vs closest plasmid homolog distance summary (.csv)
    -OUTPUT2: CDF of patristic distances in .png format
    -OUTPUT3: CDF of patristic distances in .svg format
    -OUTPUT4: CDF of patristic distances in .pdf format
"""


file1 = "Plasmid_Chr_TssC_HPC_99id_100cov.csv"
file2 = "Patristic_Distance_TssC_Matrix.tsv" 
file3 = "Co-resident_Chr_Plasmids.csv"
file4 = "Minimal_patristic_distance_CHR-PLA.csv"

OUTPUT1 = "Patristic_distances_coresident_Plasmid_Chr_TssC.csv"
OUTPUT2 = "CDF_Patristic_distances_coresident_Plasmid_Chr_TssC.png"
OUTPUT3 = "CDF_Patristic_distances_coresident_Plasmid_Chr_TssC.svg"
OUTPUT4 = "CDF_Patristic_distances_coresident_Plasmid_Chr_TssC.pdf"


df1 = pd.read_csv(file1, sep=",")
df2 = pd.read_csv(file2, sep=",", index_col = 0)

# list of chromosomal genomes
list_chr = list(set(df1[df1["Genome_Type"]=="Chr"]["Genome"].to_list()))
print("Total chromosomes: ", len(list_chr))
# list of plasmid genomes
list_plasm = list(set(df1[df1["Genome_Type"]=="Pla"]["Genome"].to_list()))
print("Total plasmids: ", len(list_plasm))

# associate each plasmid with co-resident chromosome
df3 = pd.read_csv(file3, sep=",")
dicc_coresident = {}
for i, row in df3.iterrows():
    chr = row["Chromosome"]
    plasm = row["Plasmid"]
    if chr in dicc_coresident.keys():
        dicc_coresident[chr].append(plasm)
    else:
        dicc_coresident[chr] = [plasm]

# load Chr - closest plasmid distance data
df4 = pd.read_csv(file4, sep=",")

# extract patristic distances for co-resident chromosomes and plasmids and build df for CDF graph
lista_dicc=[]
for chr_acc, plasmid_list in dicc_coresident.items():
    if chr_acc not in list_chr:
        continue
    df_chr_tssC = df1[df1["Genome"]==chr_acc]
    for i, row_chr in df_chr_tssC.iterrows():
        tssC_chr = row_chr["ID"]
        cluster_chr=row_chr["Representative"]
        for acc_plasm in plasmid_list:
            if acc_plasm not in list_plasm:
                continue
            # find all TssC in the co-resident plasmid
            df_plasm_tssC = df1[df1["Genome"] == acc_plasm]                  
            for j, row_plasm in df_plasm_tssC.iterrows():
                tssC_plasm = row_plasm["ID"]
                cluster_plasm = row_plasm["Representative"]
                distance = 0
                if cluster_chr != cluster_plasm:
                    # get distance in the matrix (df2)
                    distance = df2.at[cluster_plasm, cluster_chr]
                lista_dicc.append({
                    "Chromosomal_TssC": tssC_chr,
                    "Chr_TssC_Representative": cluster_chr,
                    "Plasmid_TssC": tssC_plasm,
                    "Plasmid_TssC_Representative": cluster_plasm,
                    "Distance": distance,
                    "Distance_Type": "Co-resident plasmid"
                })
              
            # add closest plasmid homolog
            query_match = df4[df4["Representative_query"] == cluster_chr]
            if not query_match.empty:
                match_row = query_match.iloc[0]
                lista_dicc.append({
                    "Chromosomal_TssC": tssC_chr,
                    "Chr_TssC_Representative": cluster_chr,
                    "Plasmid_TssC_Representative": match_row["Representative_target"],
                    "Distance": match_row["Minimal_distance"],
                    "Distance_Type": "Closest plasmid"
                })

df_dist = pd.DataFrame(lista_dicc)
print(df_dist)
df_dist.to_csv(OUTPUT1, sep=",", index=False)

###############################################################################################
############################# CDF GRAPHS FOR PATRISTIC DISTANCES ##############################
###############################################################################################

sns.set(style="white", color_codes=True)
a4_dims = (8, 6)
fig, ax= plt.subplots(figsize=a4_dims)
dicc_color = {"Closest plasmid" : "#cb181d", "Co-resident plasmid" : "#1380A1"}

g = sns.ecdfplot(data=df_dist, x="Distance", hue="Distance_Type", linewidth=2.5, palette=dicc_color)

plt.xlim(None, 2)
ax.set_xlabel("Patristic Distance", fontsize=16)
ax.set_ylabel("CDF", fontsize=16)
plt.setp(ax.get_legend().get_texts(), fontsize="16") 
g.legend_.set_title(None) 
g.tick_params(labelsize=14)
sns.despine()
plt.tight_layout()
plt.savefig(OUTPUT2, dpi=500);
plt.savefig(OUTPUT3, format="svg")
plt.savefig(OUTPUT4, format="pdf")
plt.close()
