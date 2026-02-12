# Plasmids drive the dissemination and diversification of Type VI Secretion Systems

This study used a comprehensive bioinformatics approach to generate a systematic, genome-wide inventory of T6SSs across bacterial genomes, evaluating their prevalence, diversity, and evolutionary relationships between plasmids and chromosomes.

The bacterial genomes were downloaded from NCBI RefSeq212 database (June 2022) and the T6SS components were detected using MacSyFinder v1.0.5 (https://github.com/gem-pasteur/macsyfinder).

This repository contains:
- The [data](data) generated in the study:
    - Co-resident chromosomal and plasmid genomes in the same bacterial hosts
    - Representative TssC proteins at 99% identity, 100% coverage used in the phylogenetic tree
    - Chromosomal and plasmid TssC phylogenetic tree file
    - Matrix of patristic distances calculated from the TssC phylogenetic tree
    - Insertion sequences (IS) detected in the proximity of plasmid T6SS gene clusters
- The relevant [Scripts](Scripts) generated and used in the study:
    - Python script to parse `hmmscan --domtblout` or `hmmsearch --domtblout` results 
    - Python script to parse `blastp -outfmt "6 std qlen slen"` results 
    - Python scripts to parse the matrix of patristic distances calculated from the TssC phylogenetic tree


For detailed information on the methods used and other supplementary material, please find a preprint version on [BioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.03.692189v1).

# Citation
If you find this work useful, please consider citing us:

María del Mar Quiñonero-Coronel, M. Pilar Garcillán-Barcia. Plasmids drive the dissemination and diversification of Type VI Secretion Systems. bioRxiv 2025.12.03.692189; doi: https://doi.org/10.64898/2025.12.03.692189 


```
@article {del Mar Qui{\~n}onero-Coronel2025.12.03.692189,
	author = {del Mar Qui{\~n}onero-Coronel, Mar{\'\i}a and Garcill{\'a}n-Barcia, M. Pilar},
	title = {{\textquotedblleft}Plasmids drive the dissemination and diversification of Type VI Secretion Systems{\textquotedblright}},
	elocation-id = {2025.12.03.692189},
	year = {2025},
	doi = {10.64898/2025.12.03.692189},
	journal = {bioRxiv}
}
```