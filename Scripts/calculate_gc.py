import os
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
import math
from tqdm import tqdm
import time
import ssl


"""
This script calculates the GC content of a list of genomes using Biopython.

Input parameters:
    -Path to the file containing the list of NCBI genome accessions (file1)

Outputs:
    -.csv file with genome accession and GC content

"""


file1 = ""
OUT_DIR = ""
FILE_NAME="Genomes_RS212_GC_"

df1 = pd.read_csv(file1)
list_genomes=df1["ID"].to_list()

BATCH_SIZE = 100


# split list if there is a large amount of genomes
if len(list_genomes)>BATCH_SIZE:
    n = math.ceil(len(list_genomes) / BATCH_SIZE)
else:
    n = len(list_genomes)

print("SEARCH: ", n)


list_dic=[]
incomplete=[]
counter=0
ssl._create_default_https_context = ssl._create_unverified_context
for i in tqdm(range (n), total=n):
    print(i)
    list_aux = list_genomes[counter: counter + BATCH_SIZE]
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=','.join(list_aux))
    recs = [rec for rec in SeqIO.parse(handle, "fasta")]
    for rec in recs:
        accession = rec.id
        sequence = rec.seq.upper()
        # obtain size
        size = len(sequence)
        g_count = sequence.count("G")
        c_count = sequence.count("C")
        gc = g_count + c_count
        gc_prop = gc*100/size
        list_dic.append({
            "Accession" : accession,
            "Size_bp" : size,
            "GC_count" : gc,
            "%GC" : gc_prop,
            "AT_count" : size - gc,
            "%AT" : 100 - gc_prop        
        })


    df_aux=pd.DataFrame(list_dic)
    df_aux.to_csv(OUT_DIR + FILE_NAME + "_" + str(i) + ".csv", sep=",", index=False)
    df_incompleto=pd.DataFrame(incomplete)
    df_incompleto.to_csv(OUT_DIR + FILE_NAME + "_incomplete_" + str(i) + ".tsv", sep="\t", index=False)

    # remove previous file
    if i > 0:
        prev_csv = OUT_DIR + FILE_NAME + "_" + str(i - 1) + ".csv"
        prev_tsv = OUT_DIR + FILE_NAME + "_incomplete_" + str(i - 1) + ".tsv"
        
        # Check if they exist before deleting to avoid errors
        if os.path.exists(prev_csv):
            os.remove(prev_csv)
        if os.path.exists(prev_tsv):
            os.remove(prev_tsv)

    counter+=BATCH_SIZE
    handle.close()
    time.sleep(1)

print("Genome accessions with GC content: ", len(list_dic))
print("Genome accessions missing GC content: ", len(incomplete))
print("FINISHED !")