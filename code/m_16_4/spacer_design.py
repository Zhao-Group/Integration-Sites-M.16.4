#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re

#Function
def revcom(seq):
    answer=''
    for i in range(0,len(seq)):       
        if seq[i]== 'A' or seq[i]== 'a':
            answer ='T'+answer
        elif seq[i]=='C' or seq[i]=='c':
            answer ='G'+answer
        elif seq[i]=='G' or seq[i]=='g':
            answer ='C'+answer
        elif seq[i]=='T' or seq[i]=='t':
            answer ='A'+answer
    return answer

#Data
df = pd.read_excel('data/m_16_4/Selected_sites.xlsx')
print(df.head())

overhang_u_left = 'AAAG'
overhang_l_left = 'TAGC'

idt_spacers = []
for i in range(np.shape(df)[0]):
    idt_spacers.append(['ID_' + str(df['ID'][i]) + '_spacer_U', overhang_u_left + df['Guide Sequence'][i] , '25nm', 'STD'])
    idt_spacers.append(['ID_' + str(df['ID'][i]) + '_spacer_L', overhang_l_left + revcom(df['Guide Sequence'][i]), '25nm', 'STD'])

idt_spacers = pd.DataFrame(idt_spacers, columns = ['Name', 'Sequence', 'Scale', 'Purification'])
idt_spacers.to_excel('data/m_16_4/m_16_4_idt_spacers'+'.xlsx', index=False)