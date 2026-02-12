import os
import pandas as pd
from tqdm import tqdm


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
ID_THRESHOLD=80
COV_THRESHOLD=80
OUTPUT1= "Results_Blastp_processed.csv"
OUTPUT2= f"Results_Blastp_{ID_THRESHOLD}id_{COV_THRESHOLD}cov.csv"
OUTPUT3= f"Results_Blastp_{ID_THRESHOLD}id_{COV_THRESHOLD}cov_besthit.csv"


def calculate_cov(row):
    "This function calculates the sequence coverage for each hit regarding the subject length"
    start = int(row["s_start"])
    end = int(row["s_end"])
    total = int(row["Subject_length"])
    cov = (end - start)*100/total
    return cov

# list to storage individual dataframes
list_df = []
directory = os.fsencode(input_dir)
for file in os.listdir(directory):
    filename=os.fsdecode(file)
    df = pd.read_csv(input_dir + filename, sep="\t", header="None")
    if df.empty:
        continue
    df.columns=["Query", "Subject", "pident", "Align_length", "Mismatch", "Gapopen", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore", "Query_length", "Subject_length"]
    # add column for genome accession
    df["Genome_Accession"]=df["Query"].apply(lambda x: x.split("|")[0])
    # calculate coverage
    tqdm.pandas(desc="my bar!")
    df["Coverage"]=df.progress_apply(calculate_cov, axis=1)

    list_df.append(df)

df_all = pd.concat([list_df])
df_all.to_csv(OUTPUT1, sep=",", index=False)
# filter by % identity and % coverage
df_filt = df_all[(df_all["pident"]>=ID_THRESHOLD) & (df_all["Coverage"]>=COV_THRESHOLD)]
df_filt.to_csv(OUTPUT2, sep=",", index=False)

##################################
# select best hit for each query #
##################################
list_aux=[]
for query, group in df_filt.groupby("Query"):
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
df_unique.to_csv(OUTPUT3, sep=",", index=False)

print("FINISHED!")



