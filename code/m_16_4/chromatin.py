#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os
import multiprocessing as mp
from Bio.Blast.Applications import NcbiblastpCommandline
import matplotlib.pyplot as plt

#Data
folder_path = 'data/m_16_4/'
df = pd.read_csv(folder_path + 'output_tr_cs.csv')
ch_input_rep1_df = pd.read_csv(folder_path + 'GSM3662076_ChIP_ClsN_input_Sis_log_rep1.txt', delimiter='\t', header=None) 
ch_input_rep2_df = pd.read_csv(folder_path + 'GSM3662077_ChIP_ClsN_input_Sis_log_rep2.txt', delimiter='\t', header=None) 
ch_IP_rep1_df = pd.read_csv(folder_path + 'GSM3662063_ChIP_ClsN_IP_Sis_log_rep1.txt', delimiter='\t', header=None) 
ch_IP_rep2_df = pd.read_csv(folder_path + 'GSM3662064_ChIP_ClsN_IP_Sis_log_rep2.txt', delimiter='\t', header=None)

#Processing
ch_df = pd.DataFrame({'Result': ((ch_IP_rep1_df[2]/ch_input_rep1_df[2]) + (ch_IP_rep2_df[2]/ch_input_rep2_df[2]))/2})
print(ch_df.head())

NUM_THREADS = mp.cpu_count()

#Functions
def read_fasta(name):
    fasta_seqs = SeqIO.parse(open(name),'fasta')
    data = []
    for fasta in fasta_seqs:
        if not ('mitochondrion' in fasta.description or 'plastid' in fasta.description or 'chloroplast' in fasta.description): #Considering only nuclear chromosomes (includes plasmid/megaplasmid), removing mitochondrial and plastid/chloroplast genomes
            data.append([fasta.id, str(fasta.seq).strip().upper()])
    
    return data

def write_fasta(name, seq_to_write):
    out_file = open(name, "w")
    for i in range(1):
        out_file.write('>' + str(1) + '\n')
        out_file.write(str(seq_to_write) + '\n')
    out_file.close()

blast_path = os.getcwd() + '/ncbi-blast-2.12.0+/bin/'
ref = folder_path + 'GCF_000189555.1_ASM18955v1_genomic.fna'
db = folder_path + 'rey15a.faa'

blastdb_cmd = '{}makeblastdb -in {} -parse_seqids -dbtype nucl -out {}'.format(blast_path, ref, db)
os.system(blastdb_cmd)

mean_chromatin = []
for i in range(len(df)): 
    seq = df['Left HR'][i][-50:] + df['PAM'][i] + df['Guide Sequence'][i] + df['Right HR'][i][0:50]
    blastout = folder_path + 'blast.tab'  # BLAST output

    if df['Strand'][i] == '+':
        write_fasta(folder_path + 'chr_temp.fa', seq)
    else:
        write_fasta(folder_path + 'chr_temp.fa', str(Seq(seq).reverse_complement()))

    cmd_blastn = NcbiblastpCommandline(cmd = blast_path + 'blastn', query = folder_path + 'chr_temp.fa', out = blastout, outfmt = 6, db = db,  num_threads=NUM_THREADS)
    stdout, stderr = cmd_blastn()

    with open(blastout, 'r') as file:
        content = file.read()

    if content:
        results = pd.read_csv(blastout, sep="\t", header=None)
        headers = ['query', 'subject',
                    'pc_identity', 'aln_length', 'mismatches', 'gaps_opened',
                    'query_start', 'query_end', 'subject_start', 'subject_end',
                    'e_value', 'bitscore']

        results.columns = headers

        print(results)
        start = results['subject_start'][0] - 1 #since we will be using python index
        stop = results['subject_end'][0] - 1 #since we will be using python index
        aln_l = results['aln_length'][0]
        if aln_l > 125:
            if stop > start:
                mean_chromatin.append(np.mean(list(ch_df[start:stop+1]['Result'])))
            else:
                mean_chromatin.append(np.mean(list(ch_df[stop:start+1]['Result'])))
        else:
            print("Alignment coverage < 90%.")
            mean_chromatin.append('-')
    else:
        print("The BLAST output file is empty.")
        mean_chromatin.append('-')

df['Mean enrichment'] = mean_chromatin
plt.plot(list(ch_df['Result']))
plt.xlabel('Position on Chromosome')
plt.ylabel('ClsN enrichment')
plt.savefig(folder_path + 'ChIP_IP_by_Input.png')
pd.DataFrame(df).to_csv('data/m_16_4/output_tr_cs_mc.csv', index = False)
