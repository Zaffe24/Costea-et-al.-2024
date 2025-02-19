## SUPPLEMENTARY FIGURE 2 j-k-i

## GSEAplot of public gene signatures

#import modules needed
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import modified_GSEAplot as plot

#import gene signatures from O'Connor et al. 2024
signatures=pd.read_table('/home/zafferan/group/MicroExonator/data_analysis/OConnor_cureated_gene_signatures.tsv',header=0,index_col=False)
gene_sets = signatures.to_dict(orient='list',)

# import list of DGEs from SCENIC analysis
df=pd.read_table('/home/zafferan/group/MicroExonator/data_analysis/stemness_markers_SCENIC.tsv',header=0,index_col=False)
df.set_index('gene',inplace=True, drop=True)

## run GSEA
results = gp.prerank(rnk=df[['avg_log2FC']],
                     gene_sets=gene_sets,
                     threads=2,
                     min_size=1,
                     max_size=10000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True # see what's going on behind the scenes
                    )

res_filt=results.res2d
res_filt.sort_values('NES',ascending=False,inplace=True)      

# save GSEA output
res_filt.to_csv('GSEA_SCENIC_stemness_markers_Oconnor.tsv',index=False,header=True,sep='\t')

# relevant signatures
terms=['T CELL QUIESCENCE_Human', 'DE NOVO PREDNISONE RESISTANCE', 'T-ALL MRD+']

# plot Suppl. Figure 2 j-k-i (one per iteration)
for t in terms:
    plot.gseaplot(term=t,
        figsize=(6,4),ofname=f"GSEAplot_{t}.png",
        legend_kws={'loc': (1.2, 0)}, **results.results[t]) 
