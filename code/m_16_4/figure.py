#Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Data
df = pd.read_excel('data/m_16_4/Data_Figure1.xlsx')
#print(df.head())

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(18, 6))

# Plot the genome length as a horizontal bar
ax.barh(0, df['Chromosome Length'][0], color='lightblue', edgecolor='black', height=0.4)

# Plot vertical lines for 'Selected Sites'
selected_sites = df['Location']
selected_colors = ['green' if site == 1 else 'black' for site in df['Selected sites']]
ax.vlines(selected_sites, -0.2, 0.2, colors=selected_colors, linestyles='dashed', label='Sites')

# Draw horizontal line for 'Conservation Score'
ax.axhline(y=-0.6, color='black')

# Plot vertical lines for 'Conservation Score'
conservation_score = df['Conservation Score']
flag = 0
for loc, score in zip(selected_sites, conservation_score):
    if flag == 0:
        ax.vlines(loc, -0.6, -0.6 + score * 0.1, colors='orange', label='Conservation Score')
        flag = 1
    else: 
        ax.vlines(loc, -0.6, -0.6 + score * 0.1, colors='orange')

for i in range(len(df) - 1):
    current_left = df.at[i, 'Left Gene']
    current_right = df.at[i, 'Right Gene']

    # Check for duplicates in subsequent rows
    for j in range(i + 1, len(df)):
        if (df.at[j, 'Left Gene'] == current_left) and (df.at[j, 'Right Gene'] == current_right):
            # If duplicates exist, set 'Left Gene Expr' and 'Right Gene Expr' to 0
            df.at[j, 'Left Gene Expr'] = 1
            df.at[j, 'Right Gene Expr'] = 1

#pd.DataFrame(df).to_csv('data/m_16_4/output_test.csv', index = False)

# Draw horizontal line for 'Gene Expr'
expr_loc = -1.2
ax.axhline(y=expr_loc, color='black')

# Plot vertical lines for 'Gene Expr'
left_gene_expr = df['Left Gene Expr'].replace('-', np.nan).astype(float)
flag = 0
for loc, expr in zip(selected_sites, left_gene_expr):
    if flag == 0:
        ax.vlines(loc - 3000, expr_loc, expr_loc + np.log(expr) * 0.05, colors='green', label='Log(Expression)')
        flag = 1
    else: 
        ax.vlines(loc - 3000, expr_loc, expr_loc + np.log(expr) * 0.05, colors='green')

right_gene_expr = df['Right Gene Expr'].replace('-', np.nan).astype(float)
for loc, expr in zip(selected_sites, right_gene_expr):
    ax.vlines(loc + 3000, expr_loc, expr_loc + np.log(expr) * 0.05, colors='green')

# Set y-axis limits to ensure labels are visible
ax.set_ylim(bottom=-1.975)

'''
OriC2: 469bp; coordinates: 1—469
OriC3: 173bp; coordinates:  1261769—1261941
OriC1:  308bp;  coordinates: 1729053–1729360
'''

ori_loc = [235, 1261855, 1729206] 
ori_labels = ['oriC2', 'oriC3', 'oriC1']

for label, loc in zip(ori_labels, ori_loc):
    ax.scatter(loc, -1.8, color='red', s=50, marker='o', zorder=5)
    ax.annotate(label, (loc, -1.97), textcoords="offset points", xytext=(16,10), ha='center', fontsize=12, zorder=10)

# Draw horizontal line for 'Mean Enrichment'
ax.axhline(y=-1.8, color='black', zorder=1)

# Plot vertical lines for 'Mean Enrichment'
mean_enrichment = df['Mean enrichment'].replace('-', np.nan).astype(float)
flag = 0
for loc, enrch in zip(selected_sites, mean_enrichment):
    if flag == 0:
        ax.vlines(loc, -1.8, -1.8 + enrch * 0.1, colors='blue', label='Mean ClsN Enrichment', zorder=2)
        flag = 1
    else:
        ax.vlines(loc, -1.8, -1.8 + enrch * 0.1, colors='blue')

# Customize the plot
ax.set_yticks([-0, -0.6, -1.2, -1.8])
ax.set_yticklabels(['M.16.4', 'Genome', 'Transcriptome', 'Epigenome'])
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')

# Show the plot
plt.tight_layout()
plt.savefig('data/m_16_4/Selected_sites.png', dpi = 400)
