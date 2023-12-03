#Modules
import numpy as np
import pandas as pd

#Data
folder_path = 'data/m_16_4/'
df = pd.read_csv(folder_path + 'output_tr_cs_mc.csv')

'''
OriC2: 469bp; coordinates: 1—469
OriC3: 173bp; coordinates:  1261769——1261941
OriC1:  308bp;  coordinates: 1729053–1729360
'''
ori_loc = [0, 1261768, 1729052] #oriC2, oriC3, oriC1
ori_size = [469, 173, 308] #oriC2, oriC3, oriC1
ori_loc_circ = ori_loc
ori_loc_circ.append(df['Chromosome Length'][0]) #incorporating circularity

dist_ori = []
closest_ori = []
for i in range(len(df)):
    dist = []
    loc = df['Location'][i]
    for j in range(len(ori_loc)): 
        if loc > ori_loc_circ[j] and loc < ori_loc_circ[j+1]: 
            dist.append(np.abs(loc - ori_loc_circ[j] - ori_size[j]))
            dist.append(np.abs(loc - ori_loc_circ[j+1]))
            break

    dist_ori.append(np.min(dist))
    if j == 0:
        if np.argmin(dist) ==  0:
            closest_ori.append('oriC2')
        else:
            closest_ori.append('oriC3')
    elif j == 1:
        if np.argmin(dist) ==  0:
            closest_ori.append('oriC3')
        else:
            closest_ori.append('oriC1')
    else: 
        if np.argmin(dist) ==  0:
            closest_ori.append('oriC1')
        else:
            closest_ori.append('oriC2')

df['Closest origin of replication'] = closest_ori
df['Distance from closest origin'] = dist_ori
pd.DataFrame(df).to_csv('data/m_16_4/output_tr_cs_mc_dO.csv', index = False)