###Supplementary Figure 2
library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(viridis)
library(readxl)
library(RColorBrewer)
library(Azimuth)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")


##Suppl Fig 2a
total.integrated <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/all_integrated/total_after_integration.rds") #load TAL1 integrated and unintegrated object
png('../../../../../../Manuscript/figures/images/Fig2/Patients_umap_unintegrated_TAL1.png',width=3500,height=2500,res=600)
DimPlot(total.integrated, group.by = "Patients", reduction = "umap.unintegrated") + ggtitle("TAL1 Patients") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Suppl Fig 2b
png('../../../../../../Manuscript/figures/images/Fig2/Patients_umap_integrated_TAL1.png',width=3500,height=2500,res=600)
DimPlot(total.integrated, group.by = "Patients", reduction = "umap") + ggtitle("TAL1 Patients") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()


##Suppl Fig 2c
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") 
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Cell_State"

prop_test <- sc_utils(obj)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Phase",
  sample_1 = "Blasts", sample_2 = "Stem_like",
  sample_identity = "Cell_State"
)
png('../../../../../../Manuscript/figures/images/Fig2/cell_cycle_permutation.png',width=4500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Stem-like - Blasts") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Phase") + theme(text =element_text(size = 15))
dev.off()

##Suppl Fig 1d-e: TAL1 alpha-beta and DP like (classification based on Pölönen et al, Nature 2024)
##alpha-beta and DP-like reference from Pölönen et al Nature 2024
TALLsubset_markers <- read_excel("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/DGEanalysis_Pölönen_Nature.xlsx")
TALLsubset_markers <- as.data.frame(TALLsubset_markers)
colnames(TALLsubset_markers)[1] <- "genes"
TAL1_DP <- subset(TALLsubset_markers, Group1 == "TAL1 DP-like")
TAL1_αβ <- subset(TALLsubset_markers, Group1 == "TAL1 αβ-like")

TAL1_DP <- TAL1_DP$genes
TAL1_αβ <- TAL1_αβ$genes


obj <- AddModuleScore(obj, features = list(TAL1_DP), name = "TAL1_DP")
obj <- AddModuleScore(obj, features = list(TAL1_αβ), name = "TAL1_αβ")

png('../../../../../../Manuscript/figures/images/Fig2/TAL1DP_Pölönen_Featureplot.png',width=3500,height=2500,res=600)
FeaturePlot(obj, "TAL1_DP1", reduction = "umap") +
  scale_color_viridis() +
  labs(
    title = "TAL1 DP-like Feature Score",
    subtitle = "Classification based on Pölönen et al, Nature 2024"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )
dev.off()

png('../../../../../../Manuscript/figures/images/Fig2/TAL1alphabeta_Pölönen_Featureplot.png',width=3500,height=2500,res=600)
FeaturePlot(obj, "TAL1_αβ1", reduction = "umap") +
  scale_color_viridis() +
  labs(
    title = "TAL1 αβ-like Feature Score",
    subtitle = "Classification based on Pölönen et al, Nature 2024"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )
dev.off()

##Suppl Fig 1f: predicted thymic celltype
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") #load TAL1 scenic object, see Figure 2
Idents(obj) <- "Cell_State" 
obj_subset <- subset(obj, subset = predicted.celltype %in% c("DN", "DP"))


prop_test <- sc_utils(obj_subset)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype",
  sample_1 = "Blasts", sample_2 = "Stem_like",
  sample_identity = "Cell_State"
)
png('../../../../../../Manuscript/figures/images/Fig2/cell_type_permutation.png',width=4500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Stem-like - Blasts") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted Thymic Celltype") + theme(text =element_text(size = 15))
dev.off()

##Suppl Fig 1g: predicted BM celltype
obj <- RunAzimuth(obj, reference = "bonemarrowref")
paired <- brewer.pal(9, "Paired")

png('../../../../../../Manuscript/figures/images/Fig2/cell_type_BM_umap.png',width=4000,height=2500,res=600)
DimPlot(obj, group.by = "predicted.celltype.l2", reduction = "umap", cols = paired) + ggtitle("Predicted Celltype using BM reference") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))  
dev.off()
##Suppl 1h
png('../../../../../../Manuscript/figures/images/Fig2/cell_type_BM_barplot.png',width=4000,height=2500,res=600)
dittoBarPlot(obj, "predicted.celltype.l2", group.by = "Cell_State", color.panel = paired)  + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + ggtitle("") + NoLegend() + xlab("Cell State")
dev.off()

##Suppl Fig 1i
obj_subset <- subset(obj, subset = predicted.celltype.l2 %in% c("CD4 Memory", "CLP"))

prop_test <- sc_utils(obj_subset)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype.l2",
  sample_1 = "Blasts", sample_2 = "Stem_like",
  sample_identity = "Cell_State"
)
png('../../../../../../Manuscript/figures/images/Fig2/BMcell_type_permutation.png',width=5500,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment Stem-like - Blasts") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Predicted BM Celltype") + theme(text =element_text(size = 15))
dev.off()




