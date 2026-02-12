import os
import pandas as pd
from tqdm import tqdm
from pathlib import Path

"""
This script processess the results of BLASTp searches (-outfmt 6) and calculate the sequence coverage for each hit.
Input parameters:
    -input_dir: Path to the directory containing .csv files from raw BLASTp search
    -ID_THRESHOLD: threshold for sequence identity
    -COV_THRESHOLD: threshold for sequence coverage
Outputs:
    -OUTPUT1: .csv file with BLASTp hits including % coverage
    -OUTPUT2: .csv file with hits using the desired thresholds
    -OUTPUT3: .csv file with best hits for each query using the desired thresholds
"""

input_dir = ""
THRESHOLD_COV = 50
VARIABLE="qlen" # "qlen" (hmmsearch: the profile is the query) or "tlen" (the profile is the target)
OUTPUT1 = "Partition_proteins_vs_PlasmidsT6SS_hmmsearch.csv"
OUTPUT2 = f"Partition_proteins_vs_PlasmidsT6SS_hmmsearch_{THRESHOLD_COV}cov.csv"
OUTPUT3 = f"Partition_proteins_vs_PlasmidsT6SS_hmmsearch_{THRESHOLD_COV}cov_besthit.csv"


def calculate_cov(row):
    "This function calculates the sequence coverage for each hit regarding the profile length"
    start = int(row["hmm_from"])
    end = int(row["hmm_to"])
    total = int(row[VARIABLE]) 
    cov = (end - start + 1)*100/total
    return cov


list_df = []

main_path = Path(input_dir)
for folder in main_path.iterdir():
    print(folder)
    folder_name = os.fsdecode(folder)
    directory=os.fsencode(folder)
    if folder.is_dir():
        for file in os.listdir(directory):
            filename=os.fsdecode(file)
            if ".domtblout" not in filename:
                continue
            with open (folder_name + "/" + filename, "r") as f:
                list_aux=[]
                rdr=iter(f)
                for line in rdr:
                    # convert each line to a dictionary
                    if line.startswith("#"):
                        continue
                    list_values=line.split()
                    dicc={
                        "Target_name" : list_values[0],
                        "Target_accession" : list_values[1],
                        "tlen" : list_values[2],
                        "Query_name": list_values[3],
                        "Query_accession" : list_values[4],
                        "qlen" : list_values[5],
                        "E-value": list_values[6],
                        "score": list_values[7], 
                        "bias": list_values[8], 
                        "#": list_values[9], 
                        "of": list_values[10],
                        "e-Evalue" : list_values[11],
                        "i-Evalue": list_values[12],
                        "score_D": list_values[13],
                        "bias_D": list_values[14],
                        "hmm_from" : list_values[15],
                        "hmm_to" : list_values[16],
                        "ali_from" : list_values[17],
                        "ali_to": list_values[18],
                        "env_from" : list_values[19],
                        "env_to" : list_values[20],
                        "acc" : list_values[21],
                        "Description_target" : " ".join(list_values[22:]),
                    }
                    list_aux.append(dicc)

            # build dataframe
            df1=pd.DataFrame(list_aux)
            if df1.empty:
                continue
            # add metadata in new columns
            df1["Genome"]=df1["Target_name"].apply(lambda x: x.split("|")[0])
            df1["CDS"]=df1["Target_name"].apply(lambda x: x.split("|")[1])
            df1["Protein_ID"]=df1["Target_name"].apply(lambda x: x.split("|")[-1])
            # calculate coverage for each hit
            tqdm.pandas(desc="my bar!")
            df1["Coverage"]=df1.progress_apply(calculate_cov, axis=1)
            list_df.append(df1)

# concat all results in one dataframe
df_all = pd.concat(list_df).reset_index(drop=True)
proteins = df_all["Target_name"].to_list()
print("Proteins identified: ", len(proteins), len(set(proteins)))
genomes = df_all["Genome"].to_list()
print("Genomes with hits: ", len(genomes), len(set(genomes)))

# filter hits by desired coverage threshold
df_filt = df_all[df_all["Coverage"]>=THRESHOLD_COV]
proteins_filt = df_filt["Target_name"].to_list()
print(f"Proteins identified >={THRESHOLD_COV}cov: ", len(proteins_filt), len(set(proteins_filt)))
genomes_filt = df_filt["Genome"].to_list()
print(f"Genomes with hits >={THRESHOLD_COV}cov: ", len(genomes_filt), len(set(genomes_filt)))
df_all.to_csv(OUTPUT1, sep="\t", index=False)
df_filt.to_csv(OUTPUT2, sep="\t", index=False)

##################################
# select best hit for each query #
##################################
list_aux=[]
for query, group in df_filt.groupby("Target_name"):
    # if there is only one heat for that protein, keep that instance
    if len(group) == 1:
        list_aux.append(group.iloc[0])
    else:
        # sort by coverage
        sorted_cov = group.sort_values(by="Coverage", ascending=False).reset_index(drop=True)
        list_aux.append(sorted_cov.iloc[0])

# obtain dataframe with only best hits
list_aux = [e.to_dict() for e in list_aux]
df_unique = pd.DataFrame(list_aux)
df_unique.to_csv(OUTPUT3, sep="\t", index=False)

print("FINISHED!")




