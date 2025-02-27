library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(org.Hs.eg.db)
library(scProportionTest)
library(ggplot2)
library(ggsignif)

##in-vitro drug treatment analysis of P1, P6 and P10 relapse
#Fig 6a)
P1.ctrl <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2024-06-05_lane16051DMSO_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P1.cyt <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2024-06-05_lane16051CYT_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)

ercc.P1.ctrl <- grep("ERCC",rownames(P1.ctrl))
ercc.P1.cyt <- grep("ERCC",rownames(P1.cyt))

P1.ctrl <- CreateSeuratObject(counts = P1.ctrl[-ercc.P1.ctrl,], project = "Control")
P1.cyt <- CreateSeuratObject(counts = P1.cyt[-ercc.P1.cyt,], project = "Cytarabine")

P1.treated <- merge(x = P1.ctrl, y = P1.cyt, add.cell.ids = c("P1_Control", "P1_Cytarabine"))

#Percentage human/mouse reads
P1.treated[["human"]] <- PercentageFeatureSet(P1.treated, pattern = "GRCh38")
P1.treated[["mouse"]] <- PercentageFeatureSet(P1.treated, pattern = "GRCm39")

#remove mouse cells + ENSEMBL gene names
P1.treated <- subset(P1.treated, subset = human > 80 )
P1.treated <- subset(P1.treated, features = rownames(P1.treated)[grepl(rownames(P1.treated),pattern = "GRCh38")])
P1.treated <- subset(P1.treated, features = rownames(P1.treated)[!grepl(rownames(P1.treated),pattern = "ENSG")])

#change gene names
gene_names <- rownames(P1.treated)
modified_gene_names <- sub("^[^-]*-", "", gene_names)
rownames(P1.treated) <- modified_gene_names
head(rownames(P1.treated))

#quality controls
P1.treated[["percent.mt"]] <- PercentageFeatureSet(P1.treated, pattern = "MT")
VlnPlot(P1.treated, features = c("percent.mt"), ncol = 1, pt.size = 0)
VlnPlot(P1.treated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
P1.treated <- subset(P1.treated, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA < 40000)

#normalization, scaling, pca
P1.treated <- NormalizeData(P1.treated)
P1.treated <- FindVariableFeatures(P1.treated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P1.treated), 10)
plot1 <- VariableFeaturePlot(P1.treated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
P1.treated <- ScaleData(P1.treated)
P1.treated <- RunPCA(P1.treated, features = VariableFeatures(object = P1.treated))

###Cell Cycle Scoring
library(readr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../../../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)


P1.treated <- JoinLayers(P1.treated)
P1.treated <- CellCycleScoring(P1.treated, 
                               g2m.features = g2m.genes, 
                               s.features = s.genes)
P1.treated <- AddModuleScore(P1.treated, features = list(cycling_histone_genes), name = "cycling_histones")
total.cell.cycle <- RunPCA(P1.treated, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")
P1.treated <- ScaleData(P1.treated, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(P1.treated))

## When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
P1.treated.cellcycle <- RunPCA(P1.treated, features = c(s.genes, g2m.genes))
DimPlot(P1.treated.cellcycle, group.by = "Phase")
FeaturePlot(P1.treated.cellcycle, features = "cycling_histones1")

##exclude noise genes + sex chromosomes (from A.Mathioudaki)
all_genes <- rownames(P1.treated@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../../../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
P1.treated <- RunPCA(object = P1.treated, verbose = FALSE, features = setdiff(VariableFeatures(object = P1.treated), exclude_genes))

ElbowPlot(P1.treated, ndims = 50)
P1.treated <- FindNeighbors(P1.treated, dims = 1:30)
P1.treated <- FindClusters(P1.treated, resolution = 0.5)
P1.treated <- RunUMAP(P1.treated, dims = 1:30)
DimPlot(P1.treated, group.by = "orig.ident")

#load stemness marker
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") #TAL1 scenic object
stemness.marker <- FindMarkers(obj, ident.1 = "Stem_like", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
stemness.marker <- subset(stemness.marker,p_val_adj < 0.05 & avg_log2FC > 0.5) 
stemness.marker <- rownames(stemness.marker)
#calculate stemness score 
P1.treated <- AddModuleScore(P1.treated, features = list(stemness.marker), name = "stemness_marker")

# Extract data from Seurat object
data <- FetchData(P1.treated, vars = c("orig.ident", "stemness_marker1"))
# Rename columns for clarity
colnames(data) <- c("Group", "Stemness_Score")
# Set the order of boxplots: Control → Cytarabine
data$Group <- factor(data$Group, levels = c("Control","Cytarabine"))
# Define comparisons for significance testing
comparisons <- list(c("Control", "Cytarabine"))
# Define custom colors for each group
custom_colors <- c("Control" = "#1f77b4", "Cytarabine" = "#ff7f0e")
# Create the boxplot
P1relinvitro.stemness <- ggplot(data, aes(x = Group, y = Stemness_Score, fill = Group)) +
  geom_boxplot(
    aes(group = Group), 
    outlier.shape = 16,   # Ensuring outliers are visible
    outlier.size = 3,     # Adjusting outlier size for clarity
    outlier.color = "black" # Outliers in black for contrast
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Explicit error bars
  geom_signif(comparisons = comparisons, 
              map_signif_level = TRUE, 
              step_increase = 0.1,  
              tip_length = 0.01,    
              color = "black") +    
  xlab("Treatment Group") + 
  ylab("Stemness Score") +
  scale_fill_manual(values = custom_colors) +  
  theme_minimal() +  
  ggtitle("P1-Relapse Stemness in-vitro") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 20, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18 , color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) +
  NoLegend()  # Remove legend for clarity

png('../../../Manuscript/figures/images/Fig6/P1rel_invitro_ctrl_cyt__barplot_stemness.png',width=4000,height=4500,res=600)
P1relinvitro.stemness
dev.off()

#Fig 6b)
P6.ctrl <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2024-06-05_lane11609DMSO_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P6.cyt <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2024-06-05_lane11609CYT_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)


ercc.P6.ctrl <- grep("ERCC",rownames(P6.ctrl))
ercc.P6.cyt <- grep("ERCC",rownames(P6.cyt))

P6.ctrl <- CreateSeuratObject(counts = P6.ctrl[-ercc.P6.ctrl,], project = "Control")
P6.cyt <- CreateSeuratObject(counts = P6.cyt[-ercc.P6.cyt,], project = "Cytarabine")

P6.treated <- merge(x = P6.ctrl, y = P6.cyt, add.cell.ids = c("P6_Control", "P6_Cytarabine"))

#Percentage human/mouse reads
P6.treated[["human"]] <- PercentageFeatureSet(P6.treated, pattern = "GRCh38")
P6.treated[["mouse"]] <- PercentageFeatureSet(P6.treated, pattern = "GRCm39")

#remove mouse cells + ENSEMBL gene names
P6.treated <- subset(P6.treated, subset = human > 80 )
P6.treated <- subset(P6.treated, features = rownames(P6.treated)[grepl(rownames(P6.treated),pattern = "GRCh38")])
P6.treated <- subset(P6.treated, features = rownames(P6.treated)[!grepl(rownames(P6.treated),pattern = "ENSG")])

#change gene names
gene_names <- rownames(P6.treated)
modified_gene_names <- sub("^[^-]*-", "", gene_names)
rownames(P6.treated) <- modified_gene_names
head(rownames(P6.treated))

#quality controls
P6.treated[["percent.mt"]] <- PercentageFeatureSet(P6.treated, pattern = "MT")
VlnPlot(P6.treated, features = c("percent.mt"), ncol = 1, pt.size = 0)
VlnPlot(P6.treated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
P6.treated <- subset(P6.treated, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA < 40000)

#normalization, scaling, pca
P6.treated <- NormalizeData(P6.treated)
P6.treated <- FindVariableFeatures(P6.treated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P6.treated), 10)
plot1 <- VariableFeaturePlot(P6.treated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
P6.treated <- ScaleData(P6.treated)
P6.treated <- RunPCA(P6.treated, features = VariableFeatures(object = P6.treated))

###Cell Cycle Scoring
library(readr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../../../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)


P6.treated <- JoinLayers(P6.treated)
P6.treated <- CellCycleScoring(P6.treated, 
                               g2m.features = g2m.genes, 
                               s.features = s.genes)
P6.treated <- AddModuleScore(P6.treated, features = list(cycling_histone_genes), name = "cycling_histones")
total.cell.cycle <- RunPCA(P6.treated, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")
P6.treated <- ScaleData(P6.treated, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(P1.treated))

## When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
P6.treated.cellcycle <- RunPCA(P6.treated, features = c(s.genes, g2m.genes))
DimPlot(P6.treated.cellcycle, group.by = "Phase")
FeaturePlot(P6.treated.cellcycle, features = "cycling_histones1")

##exclude noise genes + sex chromosomes (from A.Mathioudaki)
all_genes <- rownames(P6.treated@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../../../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
P6.treated <- RunPCA(object = P6.treated, verbose = FALSE, features = setdiff(VariableFeatures(object = P6.treated), exclude_genes))

ElbowPlot(P6.treated, ndims = 50)
P6.treated <- FindNeighbors(P6.treated, dims = 1:30)
P6.treated <- FindClusters(P6.treated, resolution = 0.5)
P6.treated <- RunUMAP(P6.treated, dims = 1:30)
DimPlot(P6.treated, group.by = "orig.ident")

#load stemness marker
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") #TAL1 scenic object
stemness.marker <- FindMarkers(obj, ident.1 = "Stem_like", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
stemness.marker <- subset(stemness.marker,p_val_adj < 0.05 & avg_log2FC > 0.5) 
stemness.marker <- rownames(stemness.marker)
#calculate stemness score 
P6.treated <- AddModuleScore(P6.treated, features = list(stemness.marker), name = "stemness_marker")

# Extract data from Seurat object
data <- FetchData(P6.treated, vars = c("orig.ident", "stemness_marker1"))
# Rename columns for clarity
colnames(data) <- c("Group", "Stemness_Score")
# Set the order of boxplots: Control → Cytarabine
data$Group <- factor(data$Group, levels = c("Control","Cytarabine"))
# Define comparisons for significance testing
comparisons <- list(c("Control", "Cytarabine"))
# Define custom colors for each group
custom_colors <- c("Control" = "#1f77b4", "Cytarabine" = "#ff7f0e")
# Create the boxplot
P6relinvitro.stemness <- ggplot(data, aes(x = Group, y = Stemness_Score, fill = Group)) +
  geom_boxplot(
    aes(group = Group), 
    outlier.shape = 16,   # Ensuring outliers are visible
    outlier.size = 3,     # Adjusting outlier size for clarity
    outlier.color = "black" # Outliers in black for contrast
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Explicit error bars
  geom_signif(comparisons = comparisons, 
              map_signif_level = TRUE, 
              step_increase = 0.1,  
              tip_length = 0.01,    
              color = "black") +    
  xlab("Treatment Group") + 
  ylab("Stemness Score") +
  scale_fill_manual(values = custom_colors) +  
  theme_minimal() +  
  ggtitle("P6-Relapse Stemness in-vitro") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 20, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18 , color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) +
  NoLegend()  # Remove legend for clarity

png('../../../Manuscript/figures/images/Fig6/P6rel_invitro_ctrl_cyt__barplot_stemness.png',width=4000,height=4500,res=600)
P6relinvitro.stemness
dev.off()

#Fig 6c)
P10rel.ctrl <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025/2025-01-21_lane1P10relapseDMSO_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P10rel.cyt <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025/2025-01-21_lane1P10relapseCYT_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)

ercc.P10.ctrl <- grep("ERCC",rownames(P10.ctrl))
ercc.P10.cyt <- grep("ERCC",rownames(P10.cyt))

P10.ctrl <- CreateSeuratObject(counts = P10.ctrl[-ercc.P10.ctrl,], project = "Control")
P10.cyt <- CreateSeuratObject(counts = P10.cyt[-ercc.P10.cyt,], project = "Cytarabine")

P10.treated <- merge(x = P10.ctrl, y = P10.cyt, add.cell.ids = c("P10_Control", "P10_Cytarabine"))

#Percentage human/mouse reads
P10.treated[["human"]] <- PercentageFeatureSet(P10.treated, pattern = "GRCh38")
P10.treated[["mouse"]] <- PercentageFeatureSet(P10.treated, pattern = "GRCm39")

#remove mouse cells + ENSEMBL gene names
P10.treated <- subset(P10.treated, subset = human > 80 )
P10.treated <- subset(P10.treated, features = rownames(P10.treated)[grepl(rownames(P10.treated),pattern = "GRCh38")])
P10.treated <- subset(P10.treated, features = rownames(P10.treated)[!grepl(rownames(P10.treated),pattern = "ENSG")])

#change gene names
gene_names <- rownames(P10.treated)
modified_gene_names <- sub("^[^-]*-", "", gene_names)
rownames(P10.treated) <- modified_gene_names
head(rownames(P10.treated))

#quality controls
P10.treated[["percent.mt"]] <- PercentageFeatureSet(P10.treated, pattern = "MT")
VlnPlot(P10.treated, features = c("percent.mt"), ncol = 1, pt.size = 0)
VlnPlot(P10.treated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
P10.treated <- subset(P10.treated, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA < 40000)

#normalization, scaling, pca
P10.treated <- NormalizeData(P10.treated)
P10.treated <- FindVariableFeatures(P10.treated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P10.treated), 10)
plot1 <- VariableFeaturePlot(P10.treated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
P10.treated <- ScaleData(P10.treated)
P10.treated <- RunPCA(P10.treated, features = VariableFeatures(object = P10.treated))

###Cell Cycle Scoring
library(readr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../../../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)


P10.treated <- JoinLayers(P10.treated)
P10.treated <- CellCycleScoring(P10.treated, 
                                g2m.features = g2m.genes, 
                                s.features = s.genes)
P10.treated <- AddModuleScore(P6.treated, features = list(cycling_histone_genes), name = "cycling_histones")
total.cell.cycle <- RunPCA(P10.treated, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")
P10.treated <- ScaleData(P10.treated, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(P1.treated))

## When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
P10.treated.cellcycle <- RunPCA(P10.treated, features = c(s.genes, g2m.genes))
DimPlot(P10.treated.cellcycle, group.by = "Phase")
FeaturePlot(P10.treated.cellcycle, features = "cycling_histones1")

##exclude noise genes + sex chromosomes (from A.Mathioudaki)
all_genes <- rownames(P10.treated@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../../../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
P10.treated <- RunPCA(object = P10.treated, verbose = FALSE, features = setdiff(VariableFeatures(object = P10.treated), exclude_genes))

ElbowPlot(P10.treated, ndims = 50)
P10.treated <- FindNeighbors(P10.treated, dims = 1:30)
P10.treated <- FindClusters(P10.treated, resolution = 0.5)
P10.treated <- RunUMAP(P10.treated, dims = 1:30)
DimPlot(P10.treated, group.by = "orig.ident")

#load stemness marker
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") #TAL1 scenic object
stemness.marker <- FindMarkers(obj, ident.1 = "Stem_like", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
stemness.marker <- subset(stemness.marker,p_val_adj < 0.05 & avg_log2FC > 0.5) 
stemness.marker <- rownames(stemness.marker)
#calculate stemness score 
P10.treated <- AddModuleScore(P10.treated, features = list(stemness.marker), name = "stemness_marker")

# Extract data from Seurat object
data <- FetchData(P10.treated, vars = c("orig.ident", "stemness_marker1"))
# Rename columns for clarity
colnames(data) <- c("Group", "Stemness_Score")
# Set the order of boxplots: Control → Cytarabine
data$Group <- factor(data$Group, levels = c("Control","Cytarabine"))
# Define comparisons for significance testing
comparisons <- list(c("Control", "Cytarabine"))
# Define custom colors for each group
custom_colors <- c("Control" = "#1f77b4", "Cytarabine" = "#ff7f0e")
# Create the boxplot
P10relinvitro.stemness <- ggplot(data, aes(x = Group, y = Stemness_Score, fill = Group)) +
  geom_boxplot(
    aes(group = Group), 
    outlier.shape = 16,   # Ensuring outliers are visible
    outlier.size = 3,     # Adjusting outlier size for clarity
    outlier.color = "black" # Outliers in black for contrast
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Explicit error bars
  geom_signif(comparisons = comparisons, 
              map_signif_level = TRUE, 
              step_increase = 0.1,  
              tip_length = 0.01,    
              color = "black") +    
  xlab("Treatment Group") + 
  ylab("Stemness Score") +
  scale_fill_manual(values = custom_colors) +  
  theme_minimal() +  
  ggtitle("P10-Relapse Stemness in-vitro") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 20, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18 , color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) +
  NoLegend()  # Remove legend for clarity

png('../../../Manuscript/figures/images/Fig6/P10rel_invitro_ctrl_cyt__barplot_stemness.png',width=4000,height=4500,res=600)
P10relinvitro.stemness
dev.off()

##in-vivo drug treatment analysis of P1 and P10 initial

#Fig 6d: P1 initial

P1ini.ctrl <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P1control_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P1ini.cyt <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P1CYT_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P1ini.mrd <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P1MRD_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)


ercc.P1ini.ctrl <- grep("ERCC",rownames(P1ini.ctrl))
ercc.P1ini.cyt <- grep("ERCC",rownames(P1ini.cyt))
ercc.P1ini.mrd <- grep("ERCC",rownames(P1ini.mrd))

P1ini.ctrl <- CreateSeuratObject(counts = P1ini.ctrl[-ercc.P1ini.ctrl,], project = "Control")
P1ini.cyt <- CreateSeuratObject(counts = P1ini.cyt[-ercc.P1ini.cyt,], project = "Cytarabine")
P1ini.mrd <- CreateSeuratObject(counts = P1ini.mrd[-ercc.P1ini.mrd,], project = "Dox-Dex-Vin")

P1ini.treated <- merge(x = P1ini.ctrl, y = c(P1ini.cyt, P1ini.mrd), add.cell.ids = c("P1-Initial_Control", "P1-Initial_Cytarabine", "P1-Initial_DoxDexVin"))

#Percentage human/mouse reads
P1ini.treated[["human"]] <- PercentageFeatureSet(P1ini.treated, pattern = "GRCh38")
P1ini.treated[["mouse"]] <- PercentageFeatureSet(P1ini.treated, pattern = "GRCm39")
head(P1ini.treated@meta.data, 5)
VlnPlot(P1ini.treated, features = c("human"))
organism.plot <- FeatureScatter(P1ini.treated, feature1 = "human", feature2 = "mouse")
organism.plot

#remove mouse cells + ENSEMBL gene names
P1ini.treated <- subset(P1ini.treated, subset = human > 80 )
P1ini.treated <- subset(P1ini.treated, features = rownames(P1ini.treated)[grepl(rownames(P1ini.treated),pattern = "GRCh38")])
P1ini.treated <- subset(P1ini.treated, features = rownames(P1ini.treated)[!grepl(rownames(P1ini.treated),pattern = "ENSG")])

#change gene names
gene_names <- rownames(P1ini.treated)
modified_gene_names <- sub("^[^-]*-", "", gene_names)
rownames(P1ini.treated) <- modified_gene_names
head(rownames(P1ini.treated))

#quality control
P1ini.treated[["percent.mt"]] <- PercentageFeatureSet(P1ini.treated, pattern = "MT")
head(P1ini.treated@meta.data, 6)
VlnPlot(P1ini.treated, features = c("percent.mt"), ncol = 1, pt.size = 0)
VlnPlot(P1ini.treated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
P1ini.treated <- subset(P1ini.treated, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA < 100000)

#normalization, dimensional reduction
P1ini.treated <- NormalizeData(P1ini.treated)
P1ini.treated <- FindVariableFeatures(P1ini.treated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P1ini.treated), 10)
plot1 <- VariableFeaturePlot(P1ini.treated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
P1ini.treated <- ScaleData(P1ini.treated)
P1ini.treated <- RunPCA(P1ini.treated, features = VariableFeatures(object = P1ini.treated))


###Cell Cycle Scoring
library(readr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)

#cell cycle scoring and scaling
P1ini.treated <- JoinLayers(P1ini.treated)
P1ini.treated <- CellCycleScoring(P1ini.treated, 
                                  g2m.features = g2m.genes, 
                                  s.features = s.genes)
P1ini.treated <- AddModuleScore(P1ini.treated, features = list(cycling_histone_genes), name = "cycling_histones")
head(P1ini.treated)
RidgePlot(P1ini.treated, features = c("S.Score", "G2M.Score", "cycling_histones1"), group.by = "orig.ident")
total.cell.cycle <- RunPCA(P1ini.treated, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

##regress out cell cycle
P1ini.treated <- ScaleData(P1ini.treated, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(P1ini.treated))

##exclude noise genes + sex chromosomes (from A.Mathioudaki)
all_genes <- rownames(P1ini.treated@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
P1ini.treated <- RunPCA(object = P1ini.treated, verbose = FALSE, features = setdiff(VariableFeatures(object = P1ini.treated), exclude_genes))

## When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
P1ini.treated.cellcycle <- RunPCA(P1ini.treated, features = c(s.genes, g2m.genes))
DimPlot(P1ini.treated.cellcycle, group.by = "Phase")
FeaturePlot(P1ini.treated.cellcycle, features = "cycling_histones1")

ElbowPlot(P1ini.treated, ndims = 50)
P1ini.treated <- FindNeighbors(P1ini.treated , dims = 1:30)
P1ini.treated  <- FindClusters(P1ini.treated , resolution = 0.5)
P1ini.treated  <- RunUMAP(P1ini.treated , dims = 1:30)
DimPlot(P1ini.treated , group.by = "orig.ident")

P1ini.treated  <- AddModuleScore(P1ini.treated , features = list(stemness.marker), name = "stemness_marker")
VlnPlot(P1ini.treated , features= "stemness_marker1", group.by="orig.ident")
saveRDS(P1ini.treated , "/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/drug_treatments/P1_Ini_in_vivo.rds")
P1ini.treated <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/drug_treatments/P1_Ini_in_vivo.rds")

# Extract data from Seurat object
data <- FetchData(P1ini.treated, vars = c("orig.ident", "stemness_marker1"))
# Rename columns for clarity
colnames(data) <- c("Group", "Stemness_Score")
# Set the order of boxplots: Control → Dox-Dex-Vin → Cytarabine
data$Group <- factor(data$Group, levels = c("Control", "Dox-Dex-Vin", "Cytarabine"))
# Define comparisons for significance testing
comparisons <- list(c("Control", "Dox-Dex-Vin"),c("Control", "Cytarabine"), c("Dox-Dex-Vin", "Cytarabine"))
# Define custom colors for each group
custom_colors <- c("Control" = "#1f77b4", "Dox-Dex-Vin" = "#2ca02c", "Cytarabine" = "#ff7f0e")

# Create the boxplot
P1ini.stemness <- ggplot(data, aes(x = Group, y = Stemness_Score, fill = Group)) +
  geom_boxplot(aes(group = Group)) +  # Boxplot grouped by category
  geom_signif(comparisons = comparisons, 
              map_signif_level = TRUE, 
              step_increase = 0.1,  # Space between brackets
              tip_length = 0.01,    # Controls bracket length
              color = "black") +    # Ensure significance annotations are black
  xlab("Treatment Group") + 
  ylab("Stemness Score") +
  scale_fill_manual(values = custom_colors) +  # Custom color scheme
  theme_minimal() +  # Removes grey background
  ggtitle("P1-Initial Stemness in-vivo") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 20, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size =  18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) +
  NoLegend()  # Remove legend to match the reference plot

png('../../../Manuscript/figures/images/Fig6/P1ini_invivo_ctrl_cyt_mrd_barplot_stemness.png',width=6000,height=4500,res=600)
print(P1ini.stemness)
dev.off()

#Fig 6e: P10 initial

#P10 Initial
P10ini.ctrl <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P10control_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P10ini.cyt <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P10CYT_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
P10ini.mrd <- Read10X("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/10x_data_TALL_drugtreated/costea/2025-02-04_lane1P10MRD_count/outs/filtered_feature_bc_matrix/", gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)

ercc.P10ini.ctrl <- grep("ERCC",rownames(P10ini.ctrl))
ercc.P10ini.cyt <- grep("ERCC",rownames(P10ini.cyt))
ercc.P10ini.mrd <- grep("ERCC",rownames(P10ini.mrd))

P10ini.ctrl <- CreateSeuratObject(counts = P10ini.ctrl[-ercc.P10ini.ctrl,], project = "Control")
P10ini.cyt <- CreateSeuratObject(counts = P10ini.cyt[-ercc.P10ini.cyt,], project = "Cytarabine")
P10ini.mrd <- CreateSeuratObject(counts = P10ini.mrd[-ercc.P10ini.mrd,], project = "Dox-Dex-Vin")

P10ini.treated <- merge(x = P10ini.ctrl, y = c(P10ini.cyt, P10ini.mrd), add.cell.ids = c("P1-Initial_Control", "P1-Initial_Cytarabine", "P1-Initial_DoxDexVin"))

#Percentage human/mouse reads
P10ini.treated[["human"]] <- PercentageFeatureSet(P10ini.treated, pattern = "GRCh38")
P10ini.treated[["mouse"]] <- PercentageFeatureSet(P10ini.treated, pattern = "GRCm39")
head(P10ini.treated@meta.data, 5)
VlnPlot(P10ini.treated, features = c("human", "mouse"), ncol = 2)
organism.plot <- FeatureScatter(P10ini.treated, feature1 = "human", feature2 = "mouse")
organism.plot

#remove mouse cells + ENSEMBL gene names
P10ini.treated <- subset(P10ini.treated, subset = human > 80 )
P10ini.treated <- subset(P10ini.treated, features = rownames(P10ini.treated)[grepl(rownames(P10ini.treated),pattern = "GRCh38")])
P10ini.treated <- subset(P10ini.treated, features = rownames(P10ini.treated)[!grepl(rownames(P10ini.treated),pattern = "ENSG")])

#change gene names
gene_names <- rownames(P10ini.treated)
modified_gene_names <- sub("^[^-]*-", "", gene_names)
rownames(P10ini.treated) <- modified_gene_names
head(rownames(P10ini.treated))

#quality control
P10ini.treated[["percent.mt"]] <- PercentageFeatureSet(P10ini.treated, pattern = "MT")
head(P10ini.treated@meta.data, 6)
VlnPlot(P10ini.treated, features = c("percent.mt"), ncol = 1, pt.size = 0)
VlnPlot(P10ini.treated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
P10ini.treated <- subset(P10ini.treated, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA < 100000)

#normalization, dimensional reduction
P10ini.treated <- NormalizeData(P10ini.treated)
P10ini.treated <- FindVariableFeatures(P10ini.treated, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P10ini.treated), 10)
plot1 <- VariableFeaturePlot(P10ini.treated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
P10ini.treated <- ScaleData(P10ini.treated)
P10ini.treated <- RunPCA(P10ini.treated, features = VariableFeatures(object = P10ini.treated))

###Cell Cycle Scoring
library(readr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)

#cell cycle scoring and scaling
P10ini.treated <- JoinLayers(P10ini.treated)
P10ini.treated <- CellCycleScoring(P10ini.treated, 
                                   g2m.features = g2m.genes, 
                                   s.features = s.genes)
P10ini.treated <- AddModuleScore(P10ini.treated, features = list(cycling_histone_genes), name = "cycling_histones")
head(P10ini.treated)
RidgePlot(P10ini.treated, features = c("S.Score", "G2M.Score", "cycling_histones1"), group.by = "orig.ident")
total.cell.cycle <- RunPCA(P10ini.treated, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")
##regress out cell cycle
P10ini.treated <- ScaleData(P10ini.treated, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(P10ini.treated))

##exclude noise genes + sex chromosomes (from A.Mathioudaki)
all_genes <- rownames(P10ini.treated@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../../../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
P10ini.treated <- RunPCA(object = P10ini.treated, verbose = FALSE, features = setdiff(VariableFeatures(object = P10ini.treated), exclude_genes))

## When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
P10ini.treated.cellcycle <- RunPCA(P10ini.treated, features = c(s.genes, g2m.genes))
DimPlot(P10ini.treated.cellcycle, group.by = "Phase")
FeaturePlot(P10ini.treated.cellcycle, features = "cycling_histones1")

ElbowPlot(P10ini.treated, ndims = 50)
P10ini.treated <- FindNeighbors(P10ini.treated , dims = 1:30)
P10ini.treated  <- FindClusters(P10ini.treated , resolution = 0.5)
P10ini.treated  <- RunUMAP(P10ini.treated , dims = 1:30)
DimPlot(P10ini.treated , group.by = "orig.ident")

P10ini.treated  <- AddModuleScore(P10ini.treated , features = list(stemness.marker), name = "stemness_marker")
VlnPlot(P10ini.treated , features= "stemness_marker1", group.by="orig.ident")
saveRDS(P10ini.treated , "/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/drug_treatments/P10_Ini_in_vivo.rds")
P10ini.treated <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/drug_treatments/P10_Ini_in_vivo.rds")

# Extract data from Seurat object
data <- FetchData(P10ini.treated, vars = c("orig.ident", "stemness_marker1"))
# Rename columns for clarity
colnames(data) <- c("Group", "Stemness_Score")
# Set the order of boxplots: Control → Dox-Dex-Vin → Cytarabine
data$Group <- factor(data$Group, levels = c("Control", "Dox-Dex-Vin", "Cytarabine"))
# Define comparisons for significance testing
comparisons <- list(c("Control", "Dox-Dex-Vin"),c("Control", "Cytarabine"), c("Dox-Dex-Vin", "Cytarabine"))
# Define custom colors for each group
custom_colors <- c("Control" = "#1f77b4", "Dox-Dex-Vin" = "#2ca02c", "Cytarabine" = "#ff7f0e")
# Create the boxplot with explicit definition of elements
P10ini.stemness <- ggplot(data, aes(x = Group, y = Stemness_Score, fill = Group)) +
  geom_boxplot(
    aes(group = Group), 
    outlier.shape = 16,   # Ensure outliers are visible
    outlier.size = 3,     # Adjusting outlier size for clarity
    outlier.color = "black" # Outliers in black for contrast
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Explicit error bars
  geom_signif(comparisons = comparisons, 
              map_signif_level = TRUE, 
              step_increase = 0.1,  
              tip_length = 0.01,    
              color = "black") +    
  xlab("Treatment Group") + 
  ylab("Stemness Score") +
  scale_fill_manual(values = custom_colors) +  
  theme_minimal() +  
  ggtitle("P10-Initial Stemness in-vivo") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 20, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size =  18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  ) +
  NoLegend()  # Remove legend for clarity

png('../../../Manuscript/figures/images/Fig6/P10ini_invivo_ctrl_cyt_mrd_barplot_stemness.png',width=6000,height=4500,res=600)
print(P10ini.stemness)
dev.off()

#Fig6f: genesetscore based on stemness marker applied to 1,300 T-ALL bulk RNA samples (Pölönen et al, Nature 2024). 
library(readxl)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(DESeq2)
library(ggpubr)


### https://github.com/HerpelinckT/geneset-modulescoring.
#This is an adaptation of Seurat's AddModuleScore() function, modified for bulk RNA sequencing data. 

AddGeneSetScore <- function(
    dds,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    name = 'Set',
    seed = 123
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  object <- counts(dds)
  object.old <- object
  object <- object %||% object.old
  
  features <- list(features)
  features.old <- features
  
  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", ")
        ) 
        warning(
          paste0("\n ",
                 paste(missing.features, collapse = ", "),
                 " dropped for calculating the geneset score."
          )
        )
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )
  
  geneset.length <- length(x = features)
  
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = object[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = geneset.length)
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = geneset.length,
    ncol = ncol(x = object)
  )
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:geneset.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  
  range01 <- lapply(
    X = features.scores.use,
    FUN = function(x) {
      range01 <- (x-min(x))/(max(x)-min(x))
    }
  )
  
  range01 <- as.data.frame(x = range01)
  rownames(x = range01) <- colnames(object)
  
  colData(dds) <- cbind(colData(dds), range01)
  
  return(dds)
}

# import metadata file of Pölönen et al, Nature 2024
meta <- read_excel("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/teacheynature2024_metadata.xlsx")
meta <- as.data.frame(meta)
meta$Day.29.morphologic.Response = factor(meta$Day.29.morphologic.Response, levels=c("M1", "M2", "M3", "Unknown"))
rownames(meta) <- meta$USI

# import expression matrix of Pölönen et al, Nature 2024
exp_mat <- data.frame(
  read.csv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/TALL_X01_counts_teacheypaper2024.tsv", 
           sep="\t"))
exp_mat$X <- sub("\\..*", "", exp_mat$X)
# converts gene IDs to gene names
gene_names <- read.csv("/g/korbel/zafferani/MicroExonator/data_analysis/gene_Id_names.tsv", sep="\t")
# not all gene Ids have names, so we cover the holes with gene Ids
gene_names$complete_Gene.name <- ifelse(gene_names$Gene.name != "", gene_names$Gene.name, gene_names$Gene.stable.ID)
# make rows unique 
ids <- gene_names %>% distinct(Gene.stable.ID,.keep_all = T)
#ids$complete_Gene.name <- make.unique(ids$complete_Gene.name)
row.names(ids) <- ids$Gene.stable.ID
#set index of expression matrix to gene names
rownames(exp_mat) <- make.names(ids[exp_mat$X, "complete_Gene.name"], unique = T)
exp_mat <- exp_mat[-c(1)]
exp_mat <- as.matrix(exp_mat)
exp_mat <- exp_mat[, rownames(meta)]

## import list of genes
genes <- read_excel("/g/korbel/Costea/Manuscript/Stemness_Marker.xlsx" )
## subset for stemness markers
## check that the thresholds correspond to the single-cell ones
stem <- subset(genes, p_val_adj < 0.05 & avg_log2FC > 0.5)
stem = as.data.frame(stem[,1])
colnames(stem) = c("complete_Gene.name")

## DESeq2
dds = DESeqDataSetFromMatrix(countData = exp_mat, colData = meta, design = ~ Day.29.morphologic.Response)
dds = AddGeneSetScore(dds, stem[,1])

comparing_groups <- list(c("M1", "M2"), c("M1", "M3"), c("M2", "M3"))
df = as.data.frame(colData(dds))
df = df[df$Day.29.morphologic.Response != "Unknown",]


png("/g/korbel/Costea/Manuscript/figures/images/Fig6/teacheybulkgenescore.png",width=5500,height=4000,res=600)
p = ggplot(data=df, aes(x=Day.29.morphologic.Response, y=Set1, fill=Day.29.morphologic.Response)) +
  geom_boxplot(
    aes(group = Day.29.morphologic.Response), 
    outlier.shape = 16,   # Ensuring outliers are visible
    outlier.size = 3,     # Adjusting outlier size for clarity
    outlier.color = "black", # Outliers in black for contrast
    width = 0.6           # Adjusting box width for better readability
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Explicit error bars
  geom_signif(comparisons = comparing_groups, 
              map_signif_level = TRUE,  # Displays stars instead of numeric p-values
              step_increase = 0.1, 
              tip_length = 0.01, 
              color = "black") +
  xlab("Day 29 morphological response") +
  ylab("Stemness score") +
  theme_minimal() +  # Removes grey background
  scale_fill_manual(values = c("M1" = "#C2C572", "M2" = "#F4A261", "M3" = "#C08497")) + # Custom colors
  NoLegend() + 
  ggtitle("Stemness-based risk profiles in 1,300 T-ALL RNA samples") +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = .5, size = 18, face = 'bold'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black")
  )
p
dev.off()
