#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
from collections import Counter
import distance
import matplotlib.pyplot as plt

#Functions
def read_fasta(name):
    fasta_seqs = SeqIO.parse(open(name),'fasta')
    data = []
    for fasta in fasta_seqs:
        if not ('mitochondrion' in fasta.description or 'plastid' in fasta.description or 'chloroplast' in fasta.description): #Considering only nuclear chromosomes, removing mitochondrial and plastid/chloroplast genomes
            data.append([fasta.id, str(fasta.seq).strip().upper()])
            
    return data
    
def pam_to_search(pam, iupac_code):
    pam_string = list(pam)
    pam_seq = []
    for i in range(len(pam_string)):
        curr_char = iupac_code[pam_string[i]].split('|')
        if i == 0:
            pam_seq = curr_char
        else:
            curr_pam_seq = []
            for j in range(len(curr_char)):
                for k in range(len(pam_seq)):
                    curr_pam_seq.append(pam_seq[k]+curr_char[j])
                    
            pam_seq = curr_pam_seq
    
    return pam_seq

def grna_search(genome, pam_l, glen, orient):
    grna_list = []
    #'for' loop for pam sequences
    for i in range(len(pam_l)):
        curr_pam = pam_l[i]
        
        #'for' loop for chromosomes
        for j in range(len(genome)):
            #top_strand
            chrom_seq = genome[j][1]
            if orient == '3prime':
                curr_grna_loc = [x - glen for x in [m.start() for m in re.finditer(curr_pam, chrom_seq)]]
            else: 
                curr_grna_loc = [m.start() for m in re.finditer(curr_pam, chrom_seq)]
                
            for k in range(len(curr_grna_loc)):
                if curr_grna_loc[k] > -1 and curr_grna_loc[k] < len(chrom_seq) - glen - len(curr_pam):
                    grna_list.append([chrom_seq[curr_grna_loc[k]:curr_grna_loc[k] + glen + len(curr_pam)], genome[j][0], curr_grna_loc[k], '+', len(chrom_seq)])
                    
            #bottom_strand
            chrom_seq = str(Seq(genome[j][1]).reverse_complement())
            if orient == '3prime':
                curr_grna_loc = [x - glen for x in [m.start() for m in re.finditer(curr_pam, chrom_seq)]]
            else: 
                curr_grna_loc = [m.start() for m in re.finditer(curr_pam, chrom_seq)]
                
            for k in range(len(curr_grna_loc)):
                if curr_grna_loc[k] > -1 and curr_grna_loc[k] < len(chrom_seq) - glen - len(curr_pam):
                    grna_list.append([chrom_seq[curr_grna_loc[k]:curr_grna_loc[k] + glen + len(curr_pam)], genome[j][0], curr_grna_loc[k], '-', len(chrom_seq)])
    
    return grna_list

iupac_code = {
  "A": "A",
  "T": "T",
  "G": "G",
  "C": "C",
  "M": "A|C",
  "R": "A|G",
  "W": "A|T",
  "S": "C|G",
  "Y": "C|T",
  "K": "G|T",
  "V": "A|C|G",
  "H": "A|C|T",
  "D": "A|G|T",
  "B": "C|G|T",
  "N": "A|C|G|T",
}

pam = 'CCN'
orient = '5prime'
glen = 40
genome = read_fasta('data/m_16_4/GCA_000022445.1_ASM2244v1_genomic.fna')
pam_library = pam_to_search(pam,iupac_code)
ambiguous_nucleotides = list(iupac_code.keys())[4:]
grna_list = grna_search(genome, pam_library, glen, orient)

#get grna occurrence table (without PAM)
if orient == '3prime':
    grna_without_pam = [item[0][0:glen] for item in grna_list]    
elif orient == '5prime':
    grna_without_pam = [item[0][len(pam):] for item in grna_list]
    
#remove ambiguous nucleotides
grna_without_pam_f = [word for word in grna_without_pam if not any(bad in word for bad in ambiguous_nucleotides)]
grna_occurrence = pd.DataFrame.from_dict(Counter(grna_without_pam_f), orient='index').reset_index()
grna_occurrence.columns = ['grna', 'frequency']

#get all guide sequences occuring in genome with duplicates removed (definition of duplicate: sequence occurring in multiple places)
complete_grna_library_wo_pam = list(grna_occurrence['grna']) 

df = pd.read_csv('data/m_16_4/output.csv')

mismatch_data = []
grna_list = list(df['Guide Sequence'])
for i in range(len(grna_list)):
    curr_d = []
    for j in range(len(complete_grna_library_wo_pam)):
        curr_d.append(distance.hamming(grna_list[i], complete_grna_library_wo_pam[j]))

    closest_mm = sorted(set(curr_d))[1]

    MM0, MM1, MM2, MM3, MM4, MM5, MM6, MM7, MM8, MM9, MM10, MM11, MM12, MM13, MM14, MM15, MM16, MM17, MM18, MM19 = [0] * 20
    if closest_mm == 0:
        MM0 = 1
    elif closest_mm == 1:
        MM1 = 1
    elif closest_mm == 2:
        MM2 = 1
    elif closest_mm == 3:
        MM3 = 1
    elif closest_mm == 4:
        MM4 = 1
    elif closest_mm == 5:
        MM5 = 1
    elif closest_mm == 6:
        MM6 = 1
    elif closest_mm == 7:
        MM7 = 1
    elif closest_mm == 8:
        MM8 = 1
    elif closest_mm == 9:
        MM9 = 1
    elif closest_mm == 10:
        MM10 = 1
    elif closest_mm == 11:
        MM11 = 1
    elif closest_mm == 12:
        MM12 = 1
    elif closest_mm == 13:
        MM13 = 1
    elif closest_mm == 14:
        MM14 = 1
    elif closest_mm == 15:
        MM15 = 1
    elif closest_mm == 16:
        MM16 = 1
    elif closest_mm == 17:
        MM17 = 1
    elif closest_mm == 18:
        MM18 = 1
    elif closest_mm == 19:
        MM19 = 1
    else:
        print('Greater or equal to 20.')
        
    if i == 0:
        mismatch_data = [[MM0, MM1, MM2, MM3, MM4, MM5, MM6, MM7, MM8, MM9, MM10, MM11, MM12, MM13, MM14, MM15, MM16, MM17, MM18, MM19]]
    else:
        mismatch_data.append([MM0, MM1, MM2, MM3, MM4, MM5, MM6, MM7, MM8, MM9, MM10, MM11, MM12, MM13, MM14, MM15, MM16, MM17, MM18, MM19])

bar_labels = ['MM0', 'MM1', 'MM2', 'MM3', 'MM4', 'MM5', 'MM6', 'MM7', 'MM8', 'MM9', 'MM10', 'MM11', 'MM12', 'MM13', 'MM14', 'MM15', 'MM16', 'MM17', 'MM18', 'MM19']
bv =  [sum(i) for i in zip(*mismatch_data)]
bar_values =  [bv[0], bv[1], bv[2], bv[3], bv[4], bv[5], bv[6], bv[7], bv[8], bv[9], bv[10], bv[11], bv[12], bv[13], bv[14], bv[15], bv[16], bv[17], bv[18], bv[19]]

fig = plt.figure(figsize = (12, 6))
bar_width = 0.4
# creating the bar plot 
plt.xlabel("Mismatches", fontsize=12)
plt.ylabel("Number of Off-targets", fontsize=12)
plt.title("Brute Force Search", fontsize=14)
bar = plt.bar(bar_labels, bar_values, bar_width, align='center')

for rect in bar:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2.0, height + 0.2, f'{height:.0f}', ha='center', va='bottom')
    
plt.ylim(top = len(grna_list))
#plt.show() 
plt.savefig('data/m_16_4/Brute Force off-target check_closest_archaea.png', dpi = 400)