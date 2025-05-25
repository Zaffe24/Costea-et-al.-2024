## python script for Supplementary Figure 4

import sys
import scanpy as sc
import pandas as pd
import os
import glob
import csv
import gzip
import pybiomart
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from decimal import Decimal
import matplotlib.cm as cm
import seaborn as sns

#### Supplementary Figure 4e: UMAPs of RPL27A (node 10-11) PSI values in 4 PDX samples

# auxiliary function
def umap_axis(axes, label=True, f=10, big=True):
    
    if big:
        #horizontal arrow
        arr_x = mpatches.FancyArrowPatch((.0, .01), (.15, 0.01),
                                    mutation_scale=15, arrowstyle='->',
                                   **{'transform':axes.transAxes})
        #vertical arrow
        arr_y=mpatches.FancyArrowPatch((0.01, 0.0), (0.01, 0.15),
                                    mutation_scale=15, arrowstyle='->',
                                   **{'transform':axes.transAxes})
        axes.add_patch(arr_x)
        axes.add_patch(arr_y)
        
        if label:
            axes.annotate("UMAP 1", (0.01, -0.02), xycoords='axes fraction', va='baseline',size=f, weight='bold')
            axes.annotate("UMAP 2", (-0.02, 0.02), xycoords='axes fraction', va='baseline',size=f, weight='bold', rotation=90)
            
### load anndata objects

p2 = sc.read_h5ad("P2_.h5ad")
p6 = sc.read_h5ad("P6_.h5ad")
p8 = sc.read_h5ad("P8_.h5ad")
p12 = sc.read_h5ad("P12_.h5ad")

###### add RPL27A_10 PSI values to each cell #####
RPL27A_10=pd.read_csv('../P2_psi/psi_RPL27A_10_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_10'],sep='\t',na_values=['NA'])
RPL27A_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_10.index]
p2.obs['RPL27A_10']=RPL27A_10.loc[obj.obs_names,'RPL27A_10']

RPL27A_10=pd.read_csv('../P6_psi/psi_RPL27A_10_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_10'],sep='\t',na_values=['NA'])
RPL27A_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_10.index]
p6.obs['RPL27A_10']=RPL27A_10.loc[obj.obs_names,'RPL27A_10']

RPL27A_10=pd.read_csv('../P8_psi/psi_RPL27A_10_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_10'],sep='\t',na_values=['NA'])
RPL27A_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_10.index]
p8.obs['RPL27A_10']=RPL27A_10.loc[obj.obs_names,'RPL27A_10']

RPL27A_10=pd.read_csv('../P12_psi/psi_RPL27A_10_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_10'],sep='\t',na_values=['NA'])
RPL27A_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_10.index]
p12.obs['RPL27A_10']=RPL27A_10.loc[obj.obs_names,'RPL27A_10']



###### add RPL27A_11 PSI values to each cell #####
RPL27A_11=pd.read_csv('../P2_psi/psi_RPL27A_11_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_11'],sep='\t',na_values=['NA'])
RPL27A_11.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_11.index]
p2.obs['RPL27A_11']=RPL27A_11.loc[obj.obs_names,'RPL27A_11']

RPL27A_11=pd.read_csv('../P6_psi/psi_RPL27A_11_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_11'],sep='\t',na_values=['NA'])
RPL27A_11.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_11.index]
p6.obs['RPL27A_11']=RPL27A_11.loc[obj.obs_names,'RPL27A_11']

RPL27A_11=pd.read_csv('../P8_psi/psi_RPL27A_11_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_11'],sep='\t',na_values=['NA'])
RPL27A_11.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_11.index]
p8.obs['RPL27A_11']=RPL27A_11.loc[obj.obs_names,'RPL27A_11']

RPL27A_11=pd.read_csv('../P12_psi/psi_RPL27A_11_sorted.tsv',header=None,index_col=0,names=['AS_id','RPL27A_11'],sep='\t',na_values=['NA'])
RPL27A_11.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in RPL27A_11.index]
p12.obs['RPL27A_11']=RPL27A_11.loc[obj.obs_names,'RPL27A_11']


###### plot UMAPS for each sample ####

event='RPL27A_10'
samples = [p2,p6,p8,p12]

for sample in samples:
    fig, ax =plt.subplots(1,1,figsize=(6,6))
    
    p=sc.pl.umap(sample,color=event,ax=ax,show=False,colorbar_loc=None,size=150,cmap='cividis',
                sort_order=True,na_color='white')
    for side in ['top', 'bottom', 'left', 'right']:
            p.spines[side].set_visible(False)
    
    p.set_xlabel(None)
    p.set_ylabel(None)
    
    umap_axis(p,label=True, big=True, f=8)
    
    ax1 = fig.add_axes([0.1, 0.3, 0.03, 0.4])  # add the top Axes
    
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=None,cmap='cividis'),cax=ax1,ticks=[0,1],ticklocation='left')
    cb.ax.set_title('Ψ')


event='RPL27A_11'
samples = [p2,p6,p8,p12]

for sample in samples:
    fig, ax =plt.subplots(1,1,figsize=(6,6))
    
    p=sc.pl.umap(sample,color=event,ax=ax,show=False,colorbar_loc=None,size=150,cmap='cividis',
                sort_order=True,na_color='white')
    for side in ['top', 'bottom', 'left', 'right']:
            p.spines[side].set_visible(False)
    
    p.set_xlabel(None)
    p.set_ylabel(None)
    
    umap_axis(p,label=True, big=True, f=8)
    
    ax1 = fig.add_axes([0.1, 0.3, 0.03, 0.4])  # add the top Axes
    
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=None,cmap='cividis'),cax=ax1,ticks=[0,1],ticklocation='left')
    cb.ax.set_title('Ψ')



######################################################################################
#####################################################################################


#### Supplementary Figure 4f: UMAPs of FOS (node 9-10) PSI values in 3 PDX samples

### load anndata objects

p2 = sc.read_h5ad("P2_.h5ad")
p8 = sc.read_h5ad("P8_.h5ad")
p10 = sc.read_h5ad("P10_.h5ad")

###### add FOS_9 PSI values to each cell #####
FOS_9=pd.read_csv('../P2_psi/psi_FOS_9_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_9'],sep='\t',na_values=['NA'])
FOS_9.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_9.index]
p2.obs['FOS_9']=FOS_9.loc[obj.obs_names,'FOS_9']

FOS_9=pd.read_csv('../P8_psi/psi_FOS_9_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_9'],sep='\t',na_values=['NA'])
FOS_9.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_9.index]
p8.obs['FOS_9']=FOS_9.loc[obj.obs_names,'FOS_9']

FOS_9=pd.read_csv('../P10_psi/psi_FOS_9_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_9'],sep='\t',na_values=['NA'])
FOS_9.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_9.index]
p10.obs['FOS_9']=FOS_9.loc[obj.obs_names,'FOS_9']

###### add FOS_10 PSI values to each cell #####

FOS_10=pd.read_csv('../P2_psi/psi_FOS_10_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_10'],sep='\t',na_values=['NA'])
FOS_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_10.index]
p2.obs['FOS_10']=FOS_10.loc[obj.obs_names,'FOS_10']

FOS_10=pd.read_csv('../P8_psi/psi_FOS_10_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_10'],sep='\t',na_values=['NA'])
FOS_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_10.index]
p8.obs['FOS_10']=FOS_10.loc[obj.obs_names,'FOS_10']

FOS_10=pd.read_csv('../P10_psi/psi_FOS_10_sorted.tsv',header=None,index_col=0,names=['AS_id','FOS_10'],sep='\t',na_values=['NA'])
FOS_10.index=[f"{i.split('_')[0]}_X{i.split('_')[2].split('.')[0]}" for i in FOS_10.index]
p10.obs['FOS_10']=FOS_10.loc[obj.obs_names,'FOS_10']

###### plot UMAPS for each sample ####

event='FOS_9'
samples = [p2,p8,p10]

for sample in samples:
    fig, ax =plt.subplots(1,1,figsize=(6,6))
    
    p=sc.pl.umap(obj,color=event,ax=ax,show=False,colorbar_loc=None,size=150,cmap='cividis',
                na_color='lightgrey')
    for side in ['top', 'bottom', 'left', 'right']:
            p.spines[side].set_visible(False)
    
    p.set_xlabel(None)
    p.set_ylabel(None)
    
    umap_axis(p,label=True, big=True, f=8)
    
    ax1 = fig.add_axes([0.1, 0.3, 0.03, 0.4])  # add the top Axes
    
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=None,cmap='cividis'),cax=ax1,ticks=[0,1],ticklocation='left')
    cb.ax.set_title('Ψ')

event='FOS_10'

for sample in samples:
    fig, ax =plt.subplots(1,1,figsize=(6,6))
    
    p=sc.pl.umap(sample,color=event,ax=ax,show=False,colorbar_loc=None,size=150,cmap='cividis',
                sort_order=True,na_color='white')
    for side in ['top', 'bottom', 'left', 'right']:
            p.spines[side].set_visible(False)
    
    p.set_xlabel(None)
    p.set_ylabel(None)
    
    umap_axis(p,label=True, big=True, f=8)
    
    ax1 = fig.add_axes([0.1, 0.3, 0.03, 0.4])  # add the top Axes
    
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=None,cmap='cividis'),cax=ax1,ticks=[0,1],ticklocation='left')
    cb.ax.set_title('Ψ')
