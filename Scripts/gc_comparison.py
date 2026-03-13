import pandas as pd
from tqdm import tqdm
import numpy as np
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt


"""
This script compares GC content of co-residen plasmids and chromosomes.

Input parameters:
    -Path to the file containing the GC metadata of plasmid genomes (file1)
    -Path to the file containing the GC metadata of chromosomal genomes (file2)
    -Path to the file containing co-resident plasmids and chromosomes (file3)

Outputs:
    -.csv file with genome accession and GC content
    -.png with boxplots of GC difference between plasmids and chromosomes for plasmids with/without T6SS or orphan islands

"""

# plasmid GC content
file1 = "RS212_plasmids_GC_T6SS_MOB.csv"
# CHR GC content
file2 = "RS212_CHR_GC.csv"
# Coresident CHR and plasmid
file3 = "Co-resident_RS212_CHR_Plasmids.csv"

OUTPUT1 = "Co-resident_CHR_Plasmid_RS212_GC_comparison.csv"
OUTPUT2 = "Co-resident_CHR_plasmid_GC_diff_T6SS_Orphan-island_MOB_boxplot.png"

df1 = pd.read_csv(file1, sep=",")
df2 = pd.read_csv(file2, sep=",")
df3 = pd.read_csv(file3, sep=",")

# merge dataframes
df_chr = df2.merge(df3, on="Chromosome")
df_merge=df_chr.merge(df1, on="Plasmid")
print(df_merge)

def calculate_dif (row):
    """This function compares the GC content of a plasmid 
    regarding its coresident chromosome"""
    gc_chr = row["Chr_%GC"]
    gc_plasmid = row["%GC"]
    diff = gc_plasmid - gc_chr
    return diff

tqdm.pandas(desc="my bar!")
df_merge["Diff_GC_Plasmid-Chr"] = df_merge.progress_apply(calculate_dif, axis=1)
df_merge["GC_diff_ABS"]=df_merge["Diff_GC_Plasmid-Chr"].apply(lambda x: abs(x))
df_merge.to_csv(OUTPUT1, sep=",", index=False)

#####################################################
########## GC COMPARISONS AND STATISTICS ############
#####################################################

# Plasmid-Chr GC difference of transmissible plasmids encoding T6SS
array_t6ss_mob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="T6SS") & (df_merge["Transmissible_plasmid"]=="Transmissible")]["Diff_GC_Plasmid-Chr"].to_list())
# Plasmid-Chr GC difference of transmissible plasmids encoding orphan islands
array_orphan_mob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="Orphan island") & (df_merge["Transmissible_plasmid"]=="Transmissible")]["Diff_GC_Plasmid-Chr"].to_list())
# Plasmid-Chr GC difference of transmissible plasmids T6SS-
array_none_mob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="Absent") & (df_merge["Transmissible_plasmid"]=="Transmissible")]["Diff_GC_Plasmid-Chr"].to_list())
# Plasmid-Chr GC difference of non-transmissible plasmids encoding T6SS
array_t6ss_nomob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="T6SS") & (df_merge["Transmissible_plasmid"]=="Non-transmissible")]["Diff_GC_Plasmid-Chr"].to_list())
# Plasmid-Chr GC difference of non-transmissible plasmids encoding orphan islands
array_orphan_nomob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="Orphan island") & (df_merge["Transmissible_plasmid"]=="Non-transmissible")]["Diff_GC_Plasmid-Chr"].to_list())
# Plasmid-Chr GC difference of non-transmissible plasmids T6SS-
array_none_nomob = np.array(df_merge[(df_merge["Putative_plasmid_T6SS"]=="Absent") & (df_merge["Transmissible_plasmid"]=="Non-transmissible")]["Diff_GC_Plasmid-Chr"].to_list())

print("\nGC content analysis: ")
print("Transmissible T6SS- Pla - Chr pairs: ", len(array_none_mob))
print("Median %GC difference Transmissible Absent: ", np.median(array_none_mob))
print("Transmissible T6SS+ Pla - Chr pairs: ", len(array_t6ss_mob))
print("Median %GC difference Transmissible T6SS: ", np.median(array_t6ss_mob))
print("Transmissible Orphan Pla - Chr pairs: ", len(array_orphan_mob))
print("Median %GC difference Transmissible Orphan: ", np.median(array_orphan_mob))
print("NON-Transmissible T6SS- Pla - Chr pairs: ", len(array_none_nomob))
print("Median %GC difference NON-Transmissible Absent: ", np.median(array_none_nomob))
print("NON-Transmissible T6SS+ Pla - Chr pairs: ", len(array_t6ss_nomob))
print("Median %GC difference NON-Transmissible T6SS: ", np.median(array_t6ss_nomob))
print("NON-Transmissible Orphan Pla - Chr pairs: ", len(array_orphan_nomob))
print("Median %GC difference NON-Transmissible Orphan: ", np.median(array_orphan_nomob))

print("Results Mann-Whitney U test NON-transmissible T6SS vs Orphan island: ", mannwhitneyu(array_t6ss_nomob, array_orphan_nomob))
print("Results Mann-Whitney U test NON-transmissible T6SS vs Absent: ", mannwhitneyu(array_t6ss_nomob, array_none_nomob))
print("Results Mann-Whitney U test NON-transmissible Orphan vs Absent: ", mannwhitneyu(array_orphan_nomob, array_none_nomob))

print("Results Mann-Whitney U test Transmissible T6SS vs Orphan island: ", mannwhitneyu(array_t6ss_mob, array_orphan_mob))
print("Results Mann-Whitney U test Transmissible T6SS vs Absent: ", mannwhitneyu(array_t6ss_mob, array_none_mob))
print("Results Mann-Whitney U test Transmissible Orphan vs Absent: ", mannwhitneyu(array_orphan_mob, array_none_mob))

print("Results Mann-Whitney U test NON-transmissible T6SS vs Transmissible T6SS: ", mannwhitneyu(array_t6ss_nomob, array_t6ss_mob))
print("Results Mann-Whitney U test NON-transmissible Orphan vs Transmissible Orphan: ", mannwhitneyu(array_orphan_nomob, array_orphan_mob))
print("Results Mann-Whitney U test NON-transmissible Absent vs Transmissible Absent: ", mannwhitneyu(array_none_nomob, array_none_mob))


############################################
###### BOXPLOT T6SS vs ORPHAN ISLAND #######
############################################

sns.set(style="white", color_codes=True)
a4_dims = (8, 6)
fig, ax= plt.subplots(figsize=a4_dims)

dicc_color = {"T6SS" : "#9372A2", "Absent" : "#BDBDBD", "Orphan island" : "#C7A000"}
boxplot=sns.boxplot(x="Transmissible_plasmid", y="Diff_GC_Plasmid-Chr", data=df_merge, hue = "Putative_plasmid_T6SS", order = ["Non-transmissible", "Transmissible"], hue_order=["Absent", "Orphan island", "T6SS"], palette=dicc_color)
ax.set_xlabel("", fontsize=18)
ax.set_ylabel("GC difference (Plasmid - Chromosome)", fontsize=18)
boxplot.tick_params(labelsize=16)

sns.despine()
plt.tight_layout()
plt.savefig(OUTPUT2, dpi=500)
