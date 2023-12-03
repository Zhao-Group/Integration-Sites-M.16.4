#Modules
import numpy as np
import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
from pymsaviz import MsaViz
import subprocess
from pyfamsa import Aligner, Sequence

#Data
folder_path = 'data/m_16_4/'
files = os.listdir(folder_path)
df = pd.read_csv(folder_path + 'output_tr.csv')

#Functions
def read_multiple_fasta(names): #Note: Code is valid if there is one chromosome in all strains
    data = {}
    for name in names:
        fasta_seqs = SeqIO.parse(open(os.path.join(folder_path, name)),'fasta')
        for fasta in fasta_seqs:
            if not ('mitochondrion' in fasta.description or 'plastid' in fasta.description or 'chloroplast' in fasta.description): #Considering only nuclear chromosomes (includes plasmid/megaplasmid), removing mitochondrial and plastid/chloroplast genomes
                if not fasta.id == 'CP001402.1':
                    data[fasta.id] = str(fasta.seq).strip().upper()
        
    return data

# Filter the .fna files
fna_files = [file for file in files if file.endswith(".fna")]
genome_dict = read_multiple_fasta(fna_files)

print(list(genome_dict.keys()))

cs_score = []
for i in range(len(df)):
    seq = df['PAM'][i] + df['Guide Sequence'][i]
    count = 0

    ofile = open(folder_path + 'temp.fa', "w")
    ofile.write(">" + 'M.16.4' + "\n" + df['Left HR'][i][-50:] + seq + df['Right HR'][i][0:50] + "\n")

    if df['Strand'][i] == '+':
        for key, value in genome_dict.items():
            if seq in value:
                loc = value.find(seq)
                ofile.write(">" + key + "\n" + value[loc-50:loc+len(seq)+50] + "\n")
                count += 1
    else:
        for key, value in genome_dict.items():
            reverse_value = str(Seq(value).reverse_complement())
            if seq in reverse_value:
                loc = reverse_value.find(seq)
                ofile.write(">" + key + "\n" + reverse_value[loc-50:loc+len(seq)+50] + "\n")
                count += 1

    ofile.close()
    cs_score.append(count)

    if count > 0:
        print(count)

        #Make MSA from fasta file
        sequences = []
        for record in SeqIO.parse(folder_path + 'temp.fa', "fasta"):
            sequences.append(Sequence(record.id.encode(), str(record.seq).encode()))

        aligner = Aligner(guide_tree="upgma")
        msa = aligner.align(sequences)

        ofile = open(folder_path + 'temp_msa.fa', "w")
        for sequence in msa:
            ofile.write(">" + sequence.id.decode() + "\n" + sequence.sequence.decode() + "\n")
        
        ofile.close()
        mv = MsaViz(folder_path + 'temp_msa.fa', color_scheme="Taylor", wrap_length=80, show_grid=True, show_consensus=True)
        mv.savefig(folder_path + 'ID_' + str(df['ID'][i]) + "_msa.png")

    #print(f"The search string appears in {count} dictionary values.")

df['Conservation Score'] = cs_score
pd.DataFrame(df).to_csv('data/m_16_4/output_tr_cs.csv', index = False)