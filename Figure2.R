library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(pals)


setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")

##preprocessing TAL1 samples
#generate sparse matrix from count tables
total.45 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v045_total.TranscriptCounts.tsv.gz")) 
total.46 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v046_total.TranscriptCounts.tsv.gz"))
total.49 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v049_total.TranscriptCounts.tsv.gz"))
total.50 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v050_total.TranscriptCounts.tsv.gz"))
total.1 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v001_total.TranscriptCounts.tsv.gz"))
total.2 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v002_total.TranscriptCounts.tsv.gz"))
total.5 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v005_total.TranscriptCounts.tsv.gz"))
total.6 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v006_total.TranscriptCounts.tsv.gz"))
total.9 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v009_total.TranscriptCounts.tsv.gz"))
total.10 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v010_total.TranscriptCounts.tsv.gz"))

total.53 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v053_total.TranscriptCounts.tsv.gz"))
total.54 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v054_total.TranscriptCounts.tsv.gz"))
total.57 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v057_total.TranscriptCounts.tsv.gz"))
total.58 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v058_total.TranscriptCounts.tsv.gz"))
total.13 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v013_total.TranscriptCounts.tsv.gz"))
total.14 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v014_total.TranscriptCounts.tsv.gz"))
total.17 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v017_total.TranscriptCounts.tsv.gz"))
total.18 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v018_total.TranscriptCounts.tsv.gz"))
total.21 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v021_total.TranscriptCounts.tsv.gz"))
total.22 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v022_total.TranscriptCounts.tsv.gz"))

total.37 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v037_total.TranscriptCounts.tsv.gz"))
total.38 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v038_total.TranscriptCounts.tsv.gz"))
total.61 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v061_total.TranscriptCounts.tsv.gz"))
total.62 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v062_total.TranscriptCounts.tsv.gz"))
total.41 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v041_total.TranscriptCounts.tsv.gz"))
total.42 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v042_total.TranscriptCounts.tsv.gz"))
total.65 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v065_total.TranscriptCounts.tsv.gz"))
total.66 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v066_total.TranscriptCounts.tsv.gz"))
total.69 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v069_total.TranscriptCounts.tsv.gz"))
total.70 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v070_total.TranscriptCounts.tsv.gz"))


total.45 <- as(total.45, "dgCMatrix")
total.46 <- as(total.46, "dgCMatrix")
total.49 <- as(total.49, "dgCMatrix")
total.50 <- as(total.50, "dgCMatrix")
total.1 <- as(total.1, "dgCMatrix")
total.2 <- as(total.2, "dgCMatrix")
total.5 <- as(total.5, "dgCMatrix")
total.6 <- as(total.6, "dgCMatrix")
total.9 <- as(total.9, "dgCMatrix")
total.10 <- as(total.10, "dgCMatrix")

total.53 <- as(total.53, "dgCMatrix")
total.54 <- as(total.54, "dgCMatrix")
total.57 <- as(total.57, "dgCMatrix")
total.58 <- as(total.58, "dgCMatrix")
total.13 <- as(total.13, "dgCMatrix")
total.14 <- as(total.14, "dgCMatrix")
total.17 <- as(total.17, "dgCMatrix")
total.18 <- as(total.18, "dgCMatrix")
total.21 <- as(total.21, "dgCMatrix")
total.22 <- as(total.22, "dgCMatrix")

total.37 <- as(total.37, "dgCMatrix")
total.38 <- as(total.38, "dgCMatrix")
total.61 <- as(total.61, "dgCMatrix")
total.62 <- as(total.62, "dgCMatrix")
total.41 <- as(total.41, "dgCMatrix")
total.42 <- as(total.42, "dgCMatrix")
total.65 <- as(total.65, "dgCMatrix")
total.66 <- as(total.66, "dgCMatrix")
total.69 <- as(total.69, "dgCMatrix")
total.70 <- as(total.70, "dgCMatrix")


ercc.total.45 <- grep("ERCC",rownames(total.45))
ercc.total.46 <- grep("ERCC",rownames(total.46))
ercc.total.49 <- grep("ERCC",rownames(total.49))
ercc.total.50 <- grep("ERCC",rownames(total.50))
ercc.total.1 <- grep("ERCC",rownames(total.1))
ercc.total.2 <- grep("ERCC",rownames(total.2))
ercc.total.5 <- grep("ERCC",rownames(total.5))
ercc.total.6 <- grep("ERCC",rownames(total.6))
ercc.total.9 <- grep("ERCC",rownames(total.9))
ercc.total.10 <- grep("ERCC",rownames(total.10))

ercc.total.53 <- grep("ERCC",rownames(total.53))
ercc.total.54 <- grep("ERCC",rownames(total.54))
ercc.total.57 <- grep("ERCC",rownames(total.57))
ercc.total.58 <- grep("ERCC",rownames(total.58))
ercc.total.13 <- grep("ERCC",rownames(total.13))
ercc.total.14 <- grep("ERCC",rownames(total.14))
ercc.total.17 <- grep("ERCC",rownames(total.17))
ercc.total.18 <- grep("ERCC",rownames(total.18))
ercc.total.21 <- grep("ERCC",rownames(total.21))
ercc.total.22 <- grep("ERCC",rownames(total.22))

ercc.total.37 <- grep("ERCC",rownames(total.37))
ercc.total.38 <- grep("ERCC",rownames(total.38))
ercc.total.61 <- grep("ERCC",rownames(total.61))
ercc.total.62 <- grep("ERCC",rownames(total.62))
ercc.total.41 <- grep("ERCC",rownames(total.41))
ercc.total.42 <- grep("ERCC",rownames(total.42))
ercc.total.65 <- grep("ERCC",rownames(total.65))
ercc.total.66 <- grep("ERCC",rownames(total.66))
ercc.total.69 <- grep("ERCC",rownames(total.69))
ercc.total.70 <- grep("ERCC",rownames(total.70))

#generate seurat objects
total.object.45 <- CreateSeuratObject(counts = total.45[-ercc.total.45,], project = "EMB-JL-v045")
total.object.46 <- CreateSeuratObject(counts = total.46[-ercc.total.46,], project = "EMB-JL-v046")
total.object.49 <- CreateSeuratObject(counts = total.49[-ercc.total.49,], project = "EMB-JL-v049")
total.object.50 <- CreateSeuratObject(counts = total.50[-ercc.total.50,], project = "EMB-JL-v050")
total.object.1 <- CreateSeuratObject(counts = total.1[-ercc.total.1,], project = "EMB-JL-v001")
total.object.2 <- CreateSeuratObject(counts = total.2[-ercc.total.2,], project = "EMB-JL-v002")
total.object.5 <- CreateSeuratObject(counts = total.5[-ercc.total.5,], project = "EMB-JL-v005")
total.object.6 <- CreateSeuratObject(counts = total.6[-ercc.total.6,], project = "EMB-JL-v006")
total.object.9 <- CreateSeuratObject(counts = total.9[-ercc.total.9,], project = "EMB-JL-v009")
total.object.10 <- CreateSeuratObject(counts = total.10[-ercc.total.10,], project = "EMB-JL-v010")

total.object.53 <- CreateSeuratObject(counts = total.53[-ercc.total.53,], project = "EMB-JL-v053")
total.object.54 <- CreateSeuratObject(counts = total.54[-ercc.total.54,], project = "EMB-JL-v054")
total.object.57 <- CreateSeuratObject(counts = total.57[-ercc.total.57,], project = "EMB-JL-v057")
total.object.58 <- CreateSeuratObject(counts = total.58[-ercc.total.58,], project = "EMB-JL-v058")
total.object.13 <- CreateSeuratObject(counts = total.13[-ercc.total.13,], project = "EMB-JL-v013")
total.object.14 <- CreateSeuratObject(counts = total.14[-ercc.total.14,], project = "EMB-JL-v014")
total.object.17 <- CreateSeuratObject(counts = total.17[-ercc.total.17,], project = "EMB-JL-v017")
total.object.18 <- CreateSeuratObject(counts = total.18[-ercc.total.18,], project = "EMB-JL-v018")
total.object.21 <- CreateSeuratObject(counts = total.21[-ercc.total.21,], project = "EMB-JL-v021")
total.object.22 <- CreateSeuratObject(counts = total.22[-ercc.total.22,], project = "EMB-JL-v022")

total.object.37 <- CreateSeuratObject(counts = total.37[-ercc.total.37,], project = "EMB-JL-v037")
total.object.38 <- CreateSeuratObject(counts = total.38[-ercc.total.38,], project = "EMB-JL-v038")
total.object.61 <- CreateSeuratObject(counts = total.61[-ercc.total.61,], project = "EMB-JL-v061")
total.object.62 <- CreateSeuratObject(counts = total.62[-ercc.total.62,], project = "EMB-JL-v062")
total.object.41 <- CreateSeuratObject(counts = total.41[-ercc.total.41,], project = "EMB-JL-v041")
total.object.42 <- CreateSeuratObject(counts = total.42[-ercc.total.42,], project = "EMB-JL-v042")
total.object.65 <- CreateSeuratObject(counts = total.65[-ercc.total.65,], project = "EMB-JL-v065")
total.object.66 <- CreateSeuratObject(counts = total.66[-ercc.total.66,], project = "EMB-JL-v066")
total.object.69 <- CreateSeuratObject(counts = total.69[-ercc.total.69,], project = "EMB-JL-v069")
total.object.70 <- CreateSeuratObject(counts = total.70[-ercc.total.70,], project = "EMB-JL-v070")


#merge seurat objects
total.object <- merge(x = total.object.45, y = c(total.object.46, total.object.57, total.object.58, total.object.1, total.object.2, total.object.5, total.object.6, total.object.9, total.object.10, total.object.53, total.object.54, total.object.49, total.object.50, total.object.13, total.object.14, total.object.17, total.object.18, total.object.21, total.object.22, total.object.37, total.object.38, total.object.61, total.object.62, total.object.41, total.object.42, total.object.65, total.object.66, total.object.69, total.object.70), add.cell.ids = c("EMB-JL-v045", "EMB-JL-v046","EMB-JL-v057", "EMB-JL-v058", "EMB-JL-v001", "EMB-JL-v002", "EMB-JL-v005", "EMB-JL-v006", "EMB-JL-v009", "EMB-JL-v010", "EMB-JL-v053", "EMB-JL-v054", "EMB-JL-v049", "EMB-JL-v050", "EMB-JL-v013", "EMB-JL-v014", "EMB-JL-v017", "EMB-JL-v018", "EMB-JL-v021", "EMB-JL-v022", "EMB-JL-v037", "EMB-JL-v038","EMB-JL-v061", "EMB-JL-v062", "EMB-JL-v041", "EMB-JL-v042", "EMB-JL-v065", "EMB-JL-v066", "EMB-JL-v069", "EMB-JL-v070" ))

#add metadata with information on disease stage and patient of origin
cluster_letters.total <- as.character((rbind(matrix('P2-Ini', 384*2,1), matrix('P6-Ini', 384*2,1), matrix('P8-Ini', 384*2,1),matrix('P10-Ini', 384*2,1),matrix('P12-Ini', 384*2,1),matrix('P2-Rel', 384*2,1),matrix('P6-Rel', 384*2,1), matrix('P8-Rel', 384*2,1), matrix('P10-Rel', 384*2,1),matrix('P12-Rel', 384*2,1),matrix('P41', 384*2,1),matrix('P52', 384*2,1), matrix('P59', 384*2,1),matrix('P63', 384*2,1),matrix('P68', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#add metadata with information on patient of origin
total.object@meta.data$Patients <- "placeholder"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P41")] <- "P41"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P52")] <- "P52"
total.object@meta.data$Patients[which(total.objects@meta.data$Patient == "P59")] <- "P59"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P63")] <- "P63"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P68")] <- "P68"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P2-Ini")] <- "P2"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P2-Rel")] <- "P2"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P6-Ini")] <- "P6"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P6-Rel")] <- "P6"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P8-Ini")] <- "P8"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P8-Rel")] <- "P8"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P10-Ini")] <- "P10"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P10-Rel")] <- "P10"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P12-Ini")] <- "P12"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P12-Rel")] <- "P12"

#Percentage human/mouse reads
total.object[["human"]] <- PercentageFeatureSet(total.object, pattern = "ENSG")
total.object[["mouse"]] <- PercentageFeatureSet(total.object, pattern = "ENSMUSG")
head(total.object@meta.data, 5)
VlnPlot(total.object, features = c("human", "mouse"), ncol = 2, group.by = "Patient")
#remove mouse cells
total <- subset(total.object, subset = human > 80 )
total <- subset(total, features = rownames(total)[grepl(rownames(total),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total$RNA) <- unique_row_names
head(total@meta.data, 6)

#Quality control and selection of cells for further analysis
total[["percent.mt"]] <- PercentageFeatureSet(total.object, pattern = "MT")
head(total@meta.data, 6)
VlnPlot(total, features = c("percent.mt"), ncol = 1, group.by = "Patient")
VlnPlot(total, features = c("nCount_RNA"), ncol = 1, group.by = "Patient")
VlnPlot(total, features = c("nFeature_RNA"), ncol = 1, group.by = "Patient")
total <- subset(total, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total <- JoinLayers(total)
total[["RNA"]] <- split(total[["RNA"]], f = total$Patient)
total <- NormalizeData(total)
total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total), 10)
plot1 <- VariableFeaturePlot(total)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
total <- ScaleData(total)
total <- RunPCA(total, features = VariableFeatures(object = total))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
###Cell Cycle Scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
head(cycling_histone_genes)
head(total@meta.data)
total <- JoinLayers(total)
total <- CellCycleScoring(total, 
                          g2m.features = g2m.genes, 
                          s.features = s.genes)
total <- AddModuleScore(total, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total[["RNA"]] <- split(total[["RNA"]], f = total$Patient)
total <- ScaleData(total, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.cell.cycle <- RunPCA(total, features = c(s.genes, g2m.genes, cycling_histone_genes))
plot1 <- DimPlot(total.cell.cycle, group.by = "Phase")
plot2 <- FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total <- RunPCA(object = total, verbose = FALSE, features = setdiff(VariableFeatures(object = total), exclude_genes))

total <- FindNeighbors(total, dims = 1:50)
total <- FindClusters(total, resolution = 0.5)
total <- RunUMAP(total, dims = 1:50)

saveRDS(total, "/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/all_integrated/total_before_integration.rds")
DimPlot(total, reduction = "umap.unintegrated")

##integration with CCA
total.integrated <- IntegrateLayers(object = total, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                    verbose = FALSE)
total.integrated[["RNA"]] <- JoinLayers(total.integrated[["RNA"]])


total.integrated <- FindNeighbors(total.integrated, reduction = "integrated.cca", dims = 1:50)
total.integrated <- FindClusters(total.integrated, resolution = 0.5)
total.integrated <- RunUMAP(total.integrated, dims = 1:50, reduction = "integrated.cca")

saveRDS(total.integrated,"/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/all_integrated/total_after_integration.rds" )


obj <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/all_integrated/total_after_integration.rds")

#add pySCENIC output to object and save it
scenic_df_wide <- read.csv("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/SCENIC_run2/output/new_aucell_mtx.tsv",
                           sep = "\t", 
                           row.names = "Cell")
colnames(scenic_df_wide) <- colnames(scenic_df_wide) %>% str_replace(pattern = fixed("..."), "")
colnames(scenic_df_wide) <- colnames(scenic_df_wide) %>% str_replace(pattern = fixed("."), "-")
all_TFs <- colnames(scenic_df_wide)
obj[["scenic"]] <- CreateAssayObject(counts = t(scenic_df_wide))

#use scenic assay for dimensional reduction and cell clustering
DefaultAssay(obj) <- "scenic"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 2.5)
obj <- RunUMAP(obj, dims = 1:50)

##Fig 2b
#UMAP patients SCENIC
png('../../../../../../Manuscript/figures/images/Fig1/Patients_umap.png',width=3500,height=2500,res=600)
DimPlot(obj, group.by = "Patients", reduction = "umap") + ggtitle("TAL1 Patients") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Fig 2c
#cell state
#add metadata with stem-like vs blasts
#Create a metadata column
obj@meta.data$Cell_State <- "Blasts"
#Change the value in "Cell_State" metadata column --> merge stem-like cell cluster
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 11)] <- "Stem_like" 
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 13)] <- "Stem_like"
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 18)] <- "Stem_like"
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 20)] <- "Stem_like"
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 21)] <- "Stem_like"
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 24)] <- "Stem_like"
obj@meta.data$Cell_State[which(obj@meta.data$seurat_clusters == 7)] <- "Stem_like"

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Cell_State"

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

png('../../../../../../Manuscript/figures/images/Fig1/cell_state_umap.png',width=3500,height=2500,res=600)
DimPlot(obj, group.by = "Cell_State", reduction = "umap", cols = c("Stem_like" = "#56B4E9", "Blasts" = "grey")) + ggtitle("Cell State") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Fig 2d-e
#Cell Cycle
png('../../../../../../Manuscript/figures/images/Fig2/cell_cycle_umap.png',width=3500,height=2500,res=600)
DimPlot(obj, reduction = "umap", group.by = "Phase", cols = c("G1" =  "#D99BBD", "G2M" = "#D5C711", "S" = "#A3A3A3"))+ ggtitle("Cell Cycle") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig2/cell_cycle_barplot.png',width=3500,height=2500,res=600)
dittoBarPlot(obj, "Phase",group.by = "Cell_State", color.panel = c("#D99BBD","#D5C711","#A3A3A3") ) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + NoLegend() + xlab("Cell State")
dev.off()

##Fig 2f: GSEA see python script

## Fig 2g-h
#cell type
thymus <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/Annotation Single Cell Atlas T cell development Science/thymus.rds") ##load reference obj generated for Fig 1h
anchor <- FindTransferAnchors(reference = thymus, query = obj,
                              dims = 1:30, reference.reduction = "pca")
obj <- MapQuery(anchorset = anchor, reference = thymus, query = obj,
                refdata = list(celltype = "Anno_level_fig1"), reference.reduction = "pca", reduction.model = "umap")
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_umap.png',width=3500,height=2500,res=600)
DimPlot(obj, group.by = "predicted.celltype", reduction = "umap", cols = c("B_naive" = "#F0E442", "CD4+Tmem" = "#56B4E9", "CD8+T" = "#009E73", "CD8+Tmem" = "#1C91D4", "CD8αα" = "#0072B2", "DN" = "#D55E00", "DP" = "#E69F00"  , "ILC3" = "#666666", "NK" = "#AD7700", "TEC(myo)" = "#007756", "TEC(neuro)" = "#CC79A7", "Treg" =  "#D5C711", "αβT(entry)" ="#005685")) + ggtitle("Predicted Celltype") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))  
dev.off()

png('../../../../../../Manuscript/figures/images/Fig1/cell_type_barplot.png',width=3500,height=2500,res=600)
dittoBarPlot(obj, "predicted.celltype", group.by = "Cell_State", color.panel = c("#F0E442","#56B4E9", "#009E73",  "#1C91D4", "#0072B2", "#D55E00",  "#E69F00"  ,  "#666666",  "#AD7700",  "#007756",  "#CC79A7",  "#D5C711", "#005685"))  + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + ggtitle("") + NoLegend() + xlab("Cell State")
dev.off()

saveRDS(obj, "/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds")
obj <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds")

##Fig 2i: differential expression analysis + volcano plotting
marker <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
marker <- subset(marker,p_val_adj < 0.05 & avg_log2FC > 0.5) 
write.xlsx(marker, row.names = TRUE, "/g/korbel/Costea/Manuscript/Nature Communications/Supplementary_Table3.xlsx")


stemness.marker.volcano <- FindMarkers(obj, ident.1 = "Stem_like", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0) #logfc threshold is set to 0 to get all data for volcano plotting in first place (for better visualization). in the plot threshold is set again to 0.5

png('../../../../../../Manuscript/figures/images/Fig3/volcanoplot_DE_TAL1.png',width=5500,height=3500,res=600)
EnhancedVolcano(stemness.marker.volcano, rownames(stemness.marker.volcano), title = "RNA Markers of Stem-like Population", subtitle = "Differential Expression: Blasts vs Stem-like Cells", titleLabSize = 23, subtitleLabSize = 17,
                selectLab = c("CD44", "ADGRE5", "CD27", "COL6A2", "AHNAK", "ITGB7", "BCL2", "BCL6", "MCL1", "NOTCH1"), boxedLabels = TRUE, legendPosition = "right",
                x ="avg_log2FC", col=c('grey', 'grey', 'grey', "red"), legendLabels = c("NS", "NS", "NS", "P-value (adj) and log2FC"),
                y ="p_val_adj", pCutoff = 0.05, FCcutoff = 0.5, pointSize = 2, labSize = 5, max.overlaps = Inf, drawConnectors = TRUE, ylab = "-Log10 P(adjusted)")
dev.off()


## Fig 2j
# SCENIC output-target genes are assigned to each TF
#Note: This is what we use for network re-construction

reg_files <- list.files("SCENIC_run2/output/regulon", 
                        pattern = ".*\\(\\+\\)\\.tsv$", 
                        full.names = T)

df_list <- list()
for (file in reg_files) {
  # the regex matches any characters except "/" that are right before a "(+).tsv" and thereby fetches the TF-names
  TF_name <- str_extract(file, "[^\\/]+(?=\\(\\+\\)\\.tsv)")
  regulon_df <- read.csv(file, sep = "\t", header = F, col.names = c("target", "count"))
  regulon_df <- mutate(regulon_df, TF = TF_name) 
  df_list[[TF_name]] <- regulon_df
}

# targene_df_raw contains all target genes for the TFs
empty_df <- data.frame(TF = character(), target = character(),  count = numeric())
targene_df_raw <- purrr::reduce(df_list, bind_rows, .init = empty_df)
# make another copy with only the target genes that were used for the activity calculation 
targene_df <- filter(targene_df_raw, count > 100) #100 out 250


#Marker genes
DefaultAssay(object = obj) <- "RNA"
DE_genes <- list()
for (p in unique(obj@meta.data$Cell_State)){
  tmp_de <- FindMarkers(object = obj, ident.1 = p)
  DE_genes[[paste0(p, "_+")]] <- rownames(tmp_de[which(tmp_de$avg_log2FC > 0.25 & tmp_de$p_val_adj < 0.000005), ])
  DE_genes[[paste0(p, "_-")]] <- rownames(tmp_de[which(tmp_de$avg_log2FC < -0.25 & tmp_de$p_val_adj < 0.000005), ])
}

pop_de <- list()
pop_de[["Blasts"]] <- unique(unlist(DE_genes[c("Blasts_+", "Blasts_-")]))
pop_de[["Stem_like"]] <- unique(unlist(DE_genes[c("Stem_like_+", "Stem_like_-")]))

intersect(pop_de[["Blasts"]], unique(targene_df$target))
sig_TFs <- list()
fisher_df1 <- list()
fisher_df2 <- list()
for (p in c("Blasts","Stem_like")){
  tmp_genes <- intersect(pop_de[[p]], unique(targene_df$target)) #Take those that are in the network
  tmp_targene_df <- targene_df[which(targene_df$target %in% tmp_genes), ] #Lets look only at all the DE genes, both conditions
  pos_genes <- intersect(tmp_genes, DE_genes[[paste0(p, "_+")]]) #cr
  neg_genes <- intersect(tmp_genes, DE_genes[[paste0(p, "_-")]])
  
  count_tfs_genes <- data.frame(TF = unique(tmp_targene_df[which(tmp_targene_df$target %in% tmp_genes), "TF"])) #how many of the condition genes these TFs target
  rownames(count_tfs_genes) <- count_tfs_genes$TF
  
  table_tfs <- data.frame() #Fisher test input
  for (t in rownames(count_tfs_genes)){
    tf_targets <- tmp_targene_df[which(tmp_targene_df$TF == t & tmp_targene_df$target %in% tmp_genes), "target"]
    count_tfs_genes[t, "pos_TF"] <- length(intersect(tf_targets, pos_genes)) + 1
    count_tfs_genes[t, "neg_TF"] <- length(intersect(tf_targets, neg_genes)) + 1
    count_tfs_genes[t, "pos_rest"] <- length(setdiff(pos_genes, tf_targets)) + 1
    count_tfs_genes[t, "neg_rest"] <- length(setdiff(neg_genes, tf_targets)) + 1
    #Contigency table
    cont_table <- data.frame()
    cont_table['pos', t] <- count_tfs_genes[t, "pos_TF"]
    cont_table['pos', 'rest'] <- count_tfs_genes[t, "pos_rest"]
    cont_table['neg', t] <- count_tfs_genes[t, "neg_TF"]
    cont_table['neg', 'rest'] <- count_tfs_genes[t, "neg_rest"]
    table_tfs[t, 'TF'] <- t
    table_tfs[t, 'Fisher_pvalue'] <- fisher.test(cont_table, alternative='two.sided', conf.int = TRUE)$p.value
    table_tfs[t, 'OR'] <- fisher.test(cont_table, alternative='two.sided', conf.int = TRUE)$estimate
  }
  table_tfs$fdr <- p.adjust(table_tfs$Fisher_pvalue, method="fdr")
  table_tfs$padj <- ifelse(table_tfs$fdr < 0.05, 
                           "<0.05", 
                           "n.s.")

  table_tfs$log2OR <- log2(table_tfs$OR)
  table_tfs <- table_tfs[order(table_tfs$log2OR), ]
  levels_plot <- rownames(table_tfs)
  write_xlsx(table_tfs, path = paste0("/g/korbel/Costea/Manuscript/tables/FisherTF_targets_", p, ".xlsx" ))
  #Plot enrichment
  
  q3 <- ggplot(table_tfs, aes(x = factor(TF, levels = levels_plot), y = log2OR, fill = padj)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          legend.position = "right") + xlab("TFs") + 
    ggtitle("Regulons of Stem-like Population") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text = element_text(size = 20)) + 
    scale_fill_manual("padj", values = c("n.s." = "#BF812D", "<0.05" = "#35978F")) 
  png('../../../../../../Manuscript/figures/images/Fig2/regulons_barplot.png',width=4500,height=2500,res=600)
  plot(q3)
  sig_TFs[[paste0(p, "_+")]] <- table_tfs[which(table_tfs$padj < 0.05 & table_tfs$log2OR > 0), "TF"]
  sig_TFs[[paste0(p, "_-")]] <- table_tfs[which(table_tfs$padj < 0.05 & table_tfs$log2OR < 0), "TF"]
  dev.off()
}

saveRDS(obj, "/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total.scenic.reduction.rds") 
