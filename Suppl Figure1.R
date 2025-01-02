###Supplementary Figure 1

library(Seurat)
library(scProportionTest)
library(RColorBrewer)
library(Azimuth)
library(ggplot2)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
total.P2 <- readRDS("Analysis/P2_no_integration/total.P2.object.rds")
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

##Suppl Fig 1a: predicted cell cycle
prop_test <- sc_utils(total.P2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Phase",
  sample_1 = "0", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_cycle_permutation_P2_0,2.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Phase") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "Phase",
  sample_1 = "0", sample_2 = "1",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_cycle_permutation_P2_1,0.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster1 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Phase") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "Phase",
  sample_1 = "1", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_cycle_permutation_P2_1,2.png',width=5000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Phase") + theme(text =element_text(size = 15))
dev.off()


##Suppl Fig 1b: predicted thymic celltype
prop_test <- sc_utils(total.P2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype",
  sample_1 = "0", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_permutation_P2_0,2.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted Thymic Celltype") + theme(text =element_text(size = 15)) + NoLegend() 
dev.off()


prop_test <- sc_utils(total.P2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype",
  sample_1 = "0", sample_2 = "1",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_permutation_P2_0,1.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster1 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted Thymic Celltype") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype",
  sample_1 = "1", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_permutation_P2_1,2.png',width=5000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted Thymic Celltype") + theme(text =element_text(size = 15))
dev.off()

##Suppl Fig 1c: batch effect Patient P2 technical replicates
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
total.P2 <- readRDS("Analysis/P2_no_integration/total.P2.object.rds")

png('../../../../../../Manuscript/figures/images/Fig1/UMAP_P2replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P2, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Suppl Fig 1d-f:predicted BM cell type 

total.P2 <- RunAzimuth(total.P2, reference = "bonemarrowref") #apply human bone marrow reference from HuBMAP consortium (https://azimuth.hubmapconsortium.org/references/human_bonemarrow/)

paired <- brewer.pal(4, "Paired")
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_BM_umap_P2.png',width=4000,height=2500,res=600)
DimPlot(total.P2, reduction = "umap", group.by = c("predicted.celltype.l2"), cols = paired)+ ggtitle("Predicted Celltype using BM reference") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))  
dev.off()
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_BM_barplot_P2.png',width=4000,height=2500,res=600)
dittoBarPlot(total.P2, "predicted.celltype.l2", group.by = "seurat_clusters", color.panel = paired) + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + ggtitle("") + NoLegend() + xlab("Seurat Clusters")
dev.off()

prop_test <- sc_utils(total.P2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype.l2",
  sample_1 = "0", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/BMcell_type_permutation_P2_0,2.png',width=3500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted BM Celltype") + theme(text =element_text(size = 15)) + NoLegend() 
dev.off()


prop_test <- sc_utils(total.P2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype.l2",
  sample_1 = "0", sample_2 = "1",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/BMcell_type_permutation_P2_0,1.png',width=3500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster1 - Cluster0") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted BM Celltype") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype.l2",
  sample_1 = "1", sample_2 = "2",
  sample_identity = "seurat_clusters"
)
png('../../../../../../Manuscript/figures/images/Fig1/BMcell_type_permutation_P2_1,2.png',width=5500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Cluster2 - Cluster1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted BM Celltype") + theme(text =element_text(size = 15))
dev.off()


##Suppl Fig 1g-i: Enrichment plots see python script


