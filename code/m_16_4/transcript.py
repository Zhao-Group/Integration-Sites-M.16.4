#Modules
import numpy as np
import pandas as pd

#Data
df = pd.read_csv('data/m_16_4/output.csv')
rna_seq = pd.read_excel('data/m_16_4/RNA-seq RJW004_Changyi.5.25.2023.xlsx')[['Synonym','Expression RJW004']]
print(rna_seq.head())

#Adding transcriptomic information
left_tr = []
right_tr = []
for i in range(np.shape(df)[0]):
    l_tr = rna_seq.loc[rna_seq['Synonym'] == df['Left Gene'][i], 'Expression RJW004']
    r_tr = rna_seq.loc[rna_seq['Synonym'] == df['Right Gene'][i], 'Expression RJW004']
    
    if not l_tr.empty:
        left_tr.append(l_tr.iloc[0])
    else:
        left_tr.append('-')
        
    if not r_tr.empty:
        right_tr.append(r_tr.iloc[0])
    else:
        right_tr.append('-')
        
df['Left Gene Expr'] = left_tr
df['Right Gene Expr'] = right_tr

pd.DataFrame(df).to_csv('data/m_16_4/output_tr.csv', index = False)