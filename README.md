# Integration-Sites-M.16.4

### THE PAPER!
This repository accompanies the work ["Characterization of chromosomal integration sites for gene expression in thermoacidophilic archaea *Sulfolobus islandicus* M.16.4"](https://www.google.com).

Command for the CRISPR-COPIES pipeline to obtain genome-wide integration sites for *Sulfolobus islandicus* M.16.4*- 
```
python code/main.py -g ../data/m_16_4/GCA_000022445.1_ASM2244v1_genomic.fna -t ../data/m_16_4/GCA_000022445.1_ASM2244v1_feature_table.txt -p CCN -o 5prime -l 40 -sl 8 --edit_dist 10 --intspace 250 -out ../data/m_16_4/output.csv --distal_end_len 20000 -hr_l 500 --protein_file ../data/m_16_4/GCA_000022445.1_ASM2244v1_protein.faa --blast_org 'Sulfolobus islandicus M.16.4' --GC_grna 25,75 --polyG_grna 6 --polyT_grna 6 --RE_grna ACCTGC,GTCGAC,CGGCCG --RE_hr ACCTGC,GTCGAC,CGGCCG
```

### Reference
<details>
<summary>If you using these scripts, please cite us:</summary>

```bibtex

```
</details>
