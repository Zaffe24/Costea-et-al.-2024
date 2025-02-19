## FIGURE 1k

## GSEAplot of stemness markers based on regulon activity

#import modules needed
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#import DGE genes
df=pd.read_table('/home/zafferan/group/MicroExonator/data_analysis/stemness_markers_SCENIC.tsv',header=0,index_col=False)
df.set_index('gene',inplace=True, drop=True)

## run GSEA on pre-ranked gene list
res_df = gp.prerank(rnk=df[['avg_log2FC']],
                     gene_sets="GO_Biological_Process_2023",
                     threads=4,
                     min_size=20,
                     max_size=10000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True # see what's going on behind the scenes
                    )
                    
res_df_filt=res_df.res2d
res_df_filt.sort_values('NES',ascending=False,inplace=True)      

# save GSEA results
res_df_filt.to_csv('GSEA_SCENIC_results.tsv',sep='\t', index=False, header=True)

## relevant GO terms
terms=['Regulation Of Cell Adhesion (GO:0030155)',
'Regulation Of I-kappaB kinase/NF-kappaB Signaling (GO:0043122)',
'Transforming Growth Factor Beta Receptor Signaling Pathway (GO:0007179)',
'Chromatin Organization (GO:0006325)',
'DNA Metabolic Process (GO:0006259)',
'Positive Regulation Of Cell Cycle Process (GO:0090068)'][::-1]

#extract relevant info
hits = [res_df.results[t]['hits'] for t in terms]
runes = [res_df.results[t]['RES'] for t in terms]
# fix colors
cmap = plt.get_cmap('seismic')
num_colors = 6
colors = cmap(np.linspace(0.9, .1, num_colors))
hex_colors = [mcolors.to_hex(color) for color in colors][::-1]

# import a customized version of the plotting function to improve image resolution
#function provided with rest of scripts
import modified_GSEAplot as plot

## plot Fig 2f
plot.gseaplot2(terms=terms, RESs=runes, hits=hits, colors=hex_colors,
                figsize=(6,6),ofname='figure_2f.png',
                legend_kws={'loc': (1.2, 0)})
