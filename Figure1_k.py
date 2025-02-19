## FIGURE 1k

## GSEAplot of DGE markers of cluster 2 from P2

#import modules needed
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#import cluster 2 markers
p2=pd.read_table('/home/zafferan/group/MicroExonator/data_analysis/P2_stemness_markers_new.tsv',header=0,index_col=False)
p2.set_index('gene',inplace=True, drop=True)

## run GSEA on pre-ranked gene list
res_p2 = gp.prerank(rnk=p2[['avg_log2FC']],
                     gene_sets="GO_Biological_Process_2023",
                     threads=4,
                     min_size=5,
                     max_size=10000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True # see what's going on behind the scenes
                    )
                    
res_p2_filt=res_p2.res2d
res_p2_filt.sort_values('NES',ascending=False,inplace=True)      

# save GSEA results
res_p2_filt.to_csv('GSEA_P2_results.tsv',sep='\t', index=False, header=True)

## relevant GO terms
terms=['Cytokine-Mediated Signaling Pathway (GO:0019221)',
'Regulation Of Cell Migration (GO:0030334)',
'Plasma Membrane Organization (GO:0007009)',
'DNA Metabolic Process (GO:0006259)',
'Nucleosome Organization (GO:0034728)',
'Positive Regulation Of Cell Cycle Process (GO:0090068)'][::-1]

#extract relevant info
hits = [res_p2.results[t]['hits'] for t in terms]
runes = [res_p2.results[t]['RES'] for t in terms]
# fix colors
cmap = plt.get_cmap('seismic')
num_colors = 6
colors = cmap(np.linspace(0.9, .1, num_colors))
hex_colors = [mcolors.to_hex(color) for color in colors][::-1]

# import a customized version of the plotting function to improve image resolution
#function provided with rest of scripts
import modified_GSEAplot as plot

## plot Fig 1k
plot.gseaplot2(terms=terms, RESs=runes, hits=hits, colors=hex_colors,
                figsize=(6,6),ofname='figure_1k.png',
                legend_kws={'loc': (1.2, 0)})
