###Supplementary Figure 3

library("scProportionTest")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

obj.subgroups <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total_allsubgroups_after_integration_44k_genes.rds") ##see figure3 for scenic object containing all subgroups
prop_test <- sc_utils(obj.subgroups)

##Suppl Fig 3a: frequency of stem-like cells per T-ALL subgroup
prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "TLX1", sample_2 = "HOXA",
  sample_identity = "Subgroup")

png('../../../../../../Manuscript/figures/images/Fig3/HOXA_TLX1_permutation.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment HOXA - TLX1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()



prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "HOXA", sample_2 = "NKX2",
  sample_identity = "Subgroup"
)
png('../../../../../../Manuscript/figures/images/Fig3/HOXA_NKX2_permutation.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment NKX2 - HOXA") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "TLX1", sample_2 = "NKX2",
  sample_identity = "Subgroup")

png('../../../../../../Manuscript/figures/images/Fig3/NKX2_TLX1_permutation.png',width=5000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment NKX2 - TLX1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15))
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "HOXA", sample_2 = "TAL1",
  sample_identity = "Subgroup")

png('../../../../../../Manuscript/figures/images/Fig3/HOXA_TAL1_permutation.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment TAL1 - HOXA") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15)) + NoLegend()
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "NKX2", sample_2 = "TAL1",
  sample_identity = "Subgroup")

png('../../../../../../Manuscript/figures/images/Fig3/TAL1_NKX2_permutation.png',width=3000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) +
  ggtitle("Enrichment TAL1 - NKX2") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15)) + NoLegend() #all data points are non significant. we add the grey color to stay consistent with previous plots. otherwise they would be coloured in red
dev.off()


prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "TLX1", sample_2 = "TAL1",
  sample_identity = "Subgroup")

png('../../../../../../Manuscript/figures/images/Fig3/TAL1_TLX1_permutation.png',width=5000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment TAL1 - TLX1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15)) 
dev.off()


##Suppl Fig 3b: frequency of stem-like cells per relapse type
prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_State",
  sample_1 = "type-1", sample_2 = "type-2",
  sample_identity = "Relapse_type"
)

png('../../../../../../Manuscript/figures/images/Fig3/relapse_type_permutation.png',width=5000,height=2500,res=600)
permutation_plot(prop_test, log2FD_threshold = 0.25, FDR_threshold = 0.05) + ggtitle("Enrichment type-2 - type-1") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cell State") + theme(text =element_text(size = 15))
dev.off()



##Suppl Fig 3c: alternative approach to stemness score: use top10 markers of TAL1 differential gene expression (DGE) analysis to analyze stemness in single cells of all subgroups

#load TAL1 object
#TAL1 object
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") 
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Cell_State"
#run DGE analysis
stemness.marker <- FindMarkers(obj, ident.1 = "Stem_like", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) 
stemness.marker <- subset(stemness.marker,p_val_adj < 0.05)

#take top10 marker and compare expression between stem like cells and blasts using cells from all subgroups
#load object containing all subgroups
obj.subgroups <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total_allsubgroups_after_integration_44k_genes.rds")
DefaultAssay(obj.subgroups) <- "RNA"
comparisons <- list(c("Stem_like", "Blasts"))
top10 <- c("SIGLEC6", "AHNAK", "PIK3R5", "KRT72", "MYO1F", "ITGB7", "LGALS1", "CD44", "KRT73", "PTPRC")
plots <- lapply(top10, function(gene) {
  VlnPlot(obj.subgroups, gene, group.by = "Cell_State", cols = c("Stem_like" = "#56B4E9", "Blasts" = "grey"), pt.size = 0) + 
    stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") + 
    NoLegend() + 
    ylim(0, 6) + 
    theme(plot.title = element_text(hjust = 0.5, face = "plain")) + 
    xlab("Cell State") + 
    theme(text = element_text(size = 15)) + 
    ggtitle(gene) 
})


combined_plot <- wrap_plots(plots, ncol = 5) + 
  plot_annotation(
    title = "Enrichment of TAL1 top10 stem-like markers in cells of all subgroups",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

png('../../../../../../Manuscript/figures/images/Fig3/DGEtop10.png',width=9000,height=5500,res=600)
print(combined_plot)
dev.off()






