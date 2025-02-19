## python script for Figure 4

import importlib
import pandas as pd
import gseapy as gp
import numpy as np
from gseapy import Msigdb
import matplotlib.pyplot as plt
import modified_GSEAplot as plot

# Figure 4d: GO over-representation analysis

# import list of genes shared by at least 2 patients undergoing AS in the stem-like population
genes=pd.read_csv('only_2_list_of_genes_shared_stemlike_patients.tsv',sep='\t',header=0)

libraries='GO_Biological_Process_2023' 

# run GO analysis
enr_All = gp.enrichr(gene_list=genes.index.to_list(),
                 gene_sets=libraries,
                 organism='human',
                 outdir=None)
                 
res_All=enr_All.results  
res_All.sort_values('Adjusted P-value',ascending=True, inplace=True)
res_All.to_csv("GO_only2_results_stemlike_patients.tsv", header=True, sep='\t')

plt.rcParams['font.size'] = 12

#filtering
new_res = res_All.loc[res_All['Adjusted P-value'] < 0.05, :].copy()
new_res.sort_values('Combined Score', inplace=True, ascending=False, ignore_index=True)

new_res=new_res.loc[:15,]
color='winter'

## Figure 4d dotplot
d=plot.dotplot(new_res,
              top_term=16,
              figsize=(4,6),
              title = 'GO analysis',
              show_ring=False, # set to False to revmove outer ring
              marker='o',
              size=10,
              cmap= color,
              ofname="Figure_4d.png"
              )


# for sashimi plots of Figure 4g and 4h check FIGURE4.sh
