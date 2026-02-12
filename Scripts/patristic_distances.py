import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt


"""
This script extracts the patristic distances in the TssC tree for the closest homolog in the same or a different 
genomic platform type (Chromosome or plasmid) and represents the Cumulative Distribution Function (CDF) of those distances. 
Inputs:
    -file1: Path to the file containing TssC clusterization 99id, 100cov representatives metadata
    -file2: Path to the file containing the patristic distance matrix
OUTPUTs:
    -OUTPUT1: Closest chromosomal homologs for plasmids distances summary
    -OUTPUT2: Closest plasmid homologs for chromosomes distances summary
    -OUTPUT3: Closest plasmid homologs for plasmids distances summary
    -OUTPUT4: Closest chromosomal homologs for chromosomes distances summary
    -OUTPUT5: CDF of patristic distances in .png format
    -OUTPUT6: CDF of patristic distances in .svg format
    -OUTPUT7: CDF of patristic distances in .pdf format
"""


file1 = "Plasmid_Chr_TssC_HPC_99id_100cov.csv"
file2 = "Patristic_Distance_TssC_Matrix.tsv" 
OUTPUT1 = "Minimal_patristic_distance_PLA-CHR.csv"
OUTPUT2 = "Minimal_patristic_distance_CHR-PLA.csv"
OUTPUT3 = "Minimal_patristic_distance_PLA-PLA.csv"
OUTPUT4 = "Minimal_patristic_distance_PLA-CHR.csv"
OUTPUT5 = "Patristic_distances_plasmid_TssC.png"
OUTPUT6 = "Patristic_distances_plasmid_TssC.svg"
OUTPUT7 = "Patristic_distances_plasmid_TssC.pdf"


df1 = pd.read_csv(file1, sep=",")
df2 = pd.read_csv(file2, sep=",", index_col = 0)

# list of chromosomal proteins
chrs = df1[df1["Genome_Type"]=="Chr"]["ID"].to_list()
print("Total chromosomal proteins: ", len(chrs))
# list of chromosomal proteins
plasm = df1[df1["Genome_Type"]=="Pla"]["ID"].to_list()
print("Total plasmid proteins: ", len(plasm))


###############################################################################################
################### MINIMAL DISTANCE DIFFERENT GENOMIC PLATFORMS (PLA - CHR) ##################
###############################################################################################

def distance_diff (protein_list, dataframe, matriz_dist, n_query, n_subject):
    """
    This function searches the minimal patristic distance for two TssC proteins encoded in different genomic
    platform type. It takes a list with protein id, a dataframe with TssC clusterization metadatada, the patristic distance matrix,
    the genomic platform query (string) and the genomic platform target (string).
    """

    list_dicc = []
    for i in protein_list:
        # retrieve TssC cluster for each protein
        repre = dataframe[dataframe["ID"]==i]["Representative"].to_list()[0]
        # if the TssC cluster contains a different genome type (CHR)
        if dataframe[dataframe["Representative"]==repre].groupby("Genome_Type").count().shape[0]>1:
            # distance = 0
            dicc ={
                "Representative_query": repre,
                "Query" : i,
                "Representative_target" : repre,
                "Genome_Type" :  n_query,
                "Minimal_distance" : 0
            }
            list_dicc.append(dicc)
        
        # if the proteins belongs to a TssC with all proteins encoded in the same genomic platform type
        else:
            # select protein in the distance matrix
            s = matriz_dist[repre]
            # sort series
            s_sorted = s.sort_values(ascending=True)
            # find the closest TssC protein encoded in a different genomic platform
            for index, value in s_sorted.items():
                if index == i:
                    continue
                # select all members contained in that TssC cluster
                members = "".join(df1[df1["Representative"]==index]["ID"].to_list())
                # if there is a member encoded in a different genomic platform, select that distance
                if n_subject in members:
                    dicc ={
                        "Representative_query": repre,
                        "Query" : i,
                        "Representative_target" : index,
                        "Genome_Type" :  n_query,
                        "Minimal_distance" : value
                    }
                    list_dicc.append(dicc)
                    break

    
    df = pd.DataFrame(list_dicc)
    return (df)

### OBTAIN MINIMAL DISTANCE PLA - CHR ###
df_pla_chr = distance_diff (plasm, df1, df2, "Pla", "Chr")
df_pla_chr["Distance_Type"] = "Plasmid - Chr"
df_pla_chr.to_csv(OUTPUT1, sep=",", index=False)
print("Search for Plasmid - Chromosome distance finished !")

df_chr_pla = distance_diff (chrs, df1, df2, "Chr", "Pla")
df_chr_pla["Distance_Type"] = "Chr - Pla"
df_chr_pla.to_csv(OUTPUT2, sep=",", index=False)
print("Search for Chromosome - Plasmid distance finished !")


###############################################################################################
########################## MINIMAL DISTANCE SAME GENOMIC PLATFORMS  ###########################
###############################################################################################
def distance_same (protein_list, dataframe, matriz_dist, n_query, n_subject):
    """
    This function searches the minimal patristic distance for two TssC proteins encoded in the same genomic
    platform type. It takes a list with protein id, a dataframe with TssC clusterization metadatada, the patristic distance matrix,
    the genomic platform query (string) and the genomic platform target (string).
    """
    list_dicc = []
    for i in protein_list:
        # retrieve TssC cluster for each protein
        repre = dataframe[dataframe["ID"]==i]["Representative"].to_list()[0]
        # if cluster contains just one sequence (TssC singletons)
        if len(dataframe[dataframe["Representative"]==repre]["ID"].to_list()) ==1:
            if dataframe[dataframe["Representative"]==repre]["Genome_Type"].to_list()[0]==n_query:
                # select protein in the distance matrix
                s = matriz_dist[repre]
                # sort distances
                s_sorted = s.sort_values(ascending=True)
                # find the closest genomic plaform
                for index, value in s_sorted.items():
                    if index == i:
                        continue
                    # select all members contained in that TssC cluster
                    members = "".join(df1[df1["Representative"]==index]["ID"].to_list())
                    # if there is a member encoded in the same genomic platform type, select that distance
                    if n_subject in members:
                        dicc ={
                            "Representative_query": repre,
                            "Query" : i,
                            "Representative_target" : index,
                            "Genome_Type" :  n_query,
                            "Minimal_distance" : value
                        }
                        list_dicc.append(dicc)
                        break            
        
        # if the cluster contains >1 TssC proteins
        else:
            aux = dataframe[dataframe["Representative"]==repre].groupby("Genome_Type").size()
            for genome_type, valor in aux.items():
                # if the cluster contains a TssC encoded in the same genomic platform type
                if genome_type==n_query:
                    # if the cluster is a singleton (1 TssC member)
                    if valor==1:
                        s = matriz_dist[repre]
                        s_sorted = s.sort_values(ascending=True)
                        for index, value in s_sorted.items():
                            if index==i:
                                continue
                            members = "".join(df1[df1["Representative"]==index]["ID"].to_list())
                            if n_subject in members:
                                dicc ={
                                    "Representative_query": repre,
                                    "Query" : i,
                                    "Representative_target" : index,
                                    "Genome_Type" :  n_query,
                                    "Minimal_distance" : value
                                }
                                list_dicc.append(dicc)
                                break            
                    # if there are >1 TssC encoded in the same genomic platform type
                    else:
                        dicc ={
                            "Representative_query": repre,
                            "Query" : i,
                            "Representative_target" : repre,
                            "Genome_Type" :  n_query,
                            "Minimal_distance" : 0
                        }
                        list_dicc.append(dicc)



    df = pd.DataFrame(list_dicc)
    return (df)

### OBTAIN MINIMAL DISTANCE PLA - PLA ###
df_pla_pla = distance_same(plasm, df1, df2, "Pla", "Pla")
df_pla_pla["Distance_Type"] = "Plasmid - Plasmid"
df_pla_pla.to_csv(OUTPUT3, sep=",", index=False)
print("Search for Plasmid - Plasmid distances finished !")

### OBTAIN MINIMAL DISTANCE CHR - CHR ###
df_chr_chr = distance_same(chrs, df1, df2, "Chr", "Chr")
df_chr_chr["Distance_Type"] = "Chr - Chr"
df_chr_chr.to_csv(OUTPUT4, sep=",", index=False)
print("Search for Chromosome - Chromosome distances finished !")

###############################################################################################
############################# CDF GRAPHS FOR PATRISTIC DISTANCES ##############################
###############################################################################################
df_all = pd.concat([df_pla_chr, df_pla_pla]).reset_index(drop=True)
df_all=df_all.sort_values("Distance_Type")

sns.set(style="white", color_codes=True)
a4_dims = (8, 6)
fig, ax= plt.subplots(figsize=a4_dims)
dicc_color = {"Plasmid - Chr" : "#525252", "Plasmid - Plasmid" : "#D95F02"}
g = sns.ecdfplot(data=df_all, x="Minimal_distance", hue="Distance_Type", linewidth=2.5, palette=dicc_color)
ax.set_xlabel("Patristic Distance", fontsize=16)
ax.set_ylabel("CDF", fontsize=16)
plt.setp(ax.get_legend().get_texts(), fontsize="16") # legend text
g.legend_.set_title(None) # remove legend title
g.tick_params(labelsize=14)
sns.despine()
plt.tight_layout()
plt.savefig(OUTPUT5, dpi=500);
plt.savefig(OUTPUT6, format="svg")
plt.savefig(OUTPUT7, format="pdf")


