library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")


##preprocessing individual patients with stem-like cells > 5%
#P6
#generate sparse matrix from count tables of patient P6 (initial and relapse)
total.49 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v049_total.TranscriptCounts.tsv.gz"))
total.50 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v050_total.TranscriptCounts.tsv.gz"))
total.57 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v057_total.TranscriptCounts.tsv.gz"))
total.58 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v058_total.TranscriptCounts.tsv.gz"))
total.49 <- as(total.49, "dgCMatrix")
total.50 <- as(total.50, "dgCMatrix")
total.57 <- as(total.57, "dgCMatrix")
total.58 <- as(total.58, "dgCMatrix")
ercc.total.49 <- grep("ERCC",rownames(total.49))
ercc.total.50 <- grep("ERCC",rownames(total.50))
ercc.total.57 <- grep("ERCC",rownames(total.57))
ercc.total.58 <- grep("ERCC",rownames(total.58))

#generate seurat objects for P6
total.object.49 <- CreateSeuratObject(counts = total.49[-ercc.total.49,], project = "EMB-JL-v049")
total.object.50 <- CreateSeuratObject(counts = total.50[-ercc.total.50,], project = "EMB-JL-v050")
total.object.57 <- CreateSeuratObject(counts = total.57[-ercc.total.57,], project = "EMB-JL-v057")
total.object.58 <- CreateSeuratObject(counts = total.58[-ercc.total.58,], project = "EMB-JL-v058")
#merge seurat objects
total.object.P6 <- merge(x = total.object.49, y = c(total.object.50, total.object.57, total.object.58), add.cell.ids = c("EMB-JL-v049", "EMB-JL-v050", "EMB-JL-v057", "EMB-JL-v058" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P6-Rel', 384*2,1),matrix('P6-Ini', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P6)
total.object.P6 <- AddMetaData(
  object = total.object.P6,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P6[["human"]] <- PercentageFeatureSet(total.object.P6, pattern = "ENSG")
total.object.P6[["mouse"]] <- PercentageFeatureSet(total.object.P6, pattern = "ENSMUSG")
head(total.object.P6@meta.data, 5)
VlnPlot(total.object.P6, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P6 <- subset(total.object.P6, subset = human > 80 )
total.P6 <- subset(total.P6, features = rownames(total.P6)[grepl(rownames(total.P6),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total.P6$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total.P6$RNA) <- unique_row_names
head(total.P6@meta.data, 6)

#Quality control and selection of cells for further analysis
total.P6[["percent.mt"]] <- PercentageFeatureSet(total.object.P6, pattern = "MT")
head(total.P6@meta.data, 6)
VlnPlot(total.P6, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P6, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P6 <- subset(total.P6, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P6 <- NormalizeData(total.P6)
total.P6 <- FindVariableFeatures(total.P6, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P6), 10)
plot1 <- VariableFeaturePlot(total.P6)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P6 <- ScaleData(total.P6)
total.P6 <- RunPCA(total.P6, features = VariableFeatures(object = total.P6))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P6 <- JoinLayers(total.P6)
total.P6 <- CellCycleScoring(total.P6, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P6 <- AddModuleScore(total.P6, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P6, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P6 <- ScaleData(total.P6, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P6))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P6.cellcycle <- RunPCA(total.P6, features = c(s.genes, g2m.genes))
DimPlot(total.P6.cellcycle, group.by = "Phase")
FeaturePlot(total.P6.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P6@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P6 <- RunPCA(object = total.P6, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P6), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P6, ndims = 50)
total.P6 <- FindNeighbors(total.P6, dims = 1:30)
total.P6 <- FindClusters(total.P6, resolution = 0.5)
total.P6 <- RunUMAP(total.P6, dims = 1:30)

#P8
#generate sparse matrix from count tables of patient P8 (initial and relapse)
total.1 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v001_total.TranscriptCounts.tsv.gz"))
total.2 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v002_total.TranscriptCounts.tsv.gz"))
total.13 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v013_total.TranscriptCounts.tsv.gz"))
total.14 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v014_total.TranscriptCounts.tsv.gz"))
total.1 <- as(total.1, "dgCMatrix")
total.2 <- as(total.2, "dgCMatrix")
total.13 <- as(total.13, "dgCMatrix")
total.14 <- as(total.14, "dgCMatrix")
ercc.total.1 <- grep("ERCC",rownames(total.1))
ercc.total.2 <- grep("ERCC",rownames(total.2))
ercc.total.13 <- grep("ERCC",rownames(total.13))
ercc.total.14 <- grep("ERCC",rownames(total.14))


#generate seurat objects for P8
total.object.1 <- CreateSeuratObject(counts = total.1[-ercc.total.1,], project = "EMB-JL-v001")
total.object.2 <- CreateSeuratObject(counts = total.2[-ercc.total.2,], project = "EMB-JL-v002")
total.object.13 <- CreateSeuratObject(counts = total.13[-ercc.total.13,], project = "EMB-JL-v013")
total.object.14 <- CreateSeuratObject(counts = total.14[-ercc.total.14,], project = "EMB-JL-v014")
#merge seurat objects
total.object.P8 <- merge(x = total.object.1, y = c(total.object.2, total.object.13, total.object.14), add.cell.ids = c("EMB-JL-v001", "EMB-JL-v002", "EMB-JL-v013", "EMB-JL-v014" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P8-Ini', 384*2,1),matrix('P8-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P8)
total.object.P8 <- AddMetaData(
  object = total.object.P8,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P8[["human"]] <- PercentageFeatureSet(total.object.P8, pattern = "ENSG")
total.object.P8[["mouse"]] <- PercentageFeatureSet(total.object.P8, pattern = "ENSMUSG")
head(total.object.P8@meta.data, 5)
VlnPlot(total.object.P8, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P8 <- subset(total.object.P8, subset = human > 80 )
total.P8 <- subset(total.P8, features = rownames(total.P8)[grepl(rownames(total.P8),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total.P8$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total.P8$RNA) <- unique_row_names
head(total.P8@meta.data, 6)

#Quality control and selection of cells for further analysis
total.P8[["percent.mt"]] <- PercentageFeatureSet(total.object.P8, pattern = "MT")
head(total.P8@meta.data, 6)
VlnPlot(total.P8, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P8, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P8 <- subset(total.P8, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P8 <- NormalizeData(total.P8)
total.P8 <- FindVariableFeatures(total.P8, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P8), 10)
plot1 <- VariableFeaturePlot(total.P8)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P8 <- ScaleData(total.P8)
total.P8 <- RunPCA(total.P8, features = VariableFeatures(object = total.P8))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P8 <- JoinLayers(total.P8)
total.P8 <- CellCycleScoring(total.P8, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P8 <- AddModuleScore(total.P8, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P8, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P8 <- ScaleData(total.P8, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P8))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P8.cellcycle <- RunPCA(total.P8, features = c(s.genes, g2m.genes))
DimPlot(total.P8.cellcycle, group.by = "Phase")
FeaturePlot(total.P8.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P8@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P8 <- RunPCA(object = total.P8, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P8), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P8, ndims = 50)
total.P8 <- FindNeighbors(total.P8, dims = 1:30)
total.P8 <- FindClusters(total.P8, resolution = 0.5)
total.P8 <- RunUMAP(total.P8, dims = 1:30)


#P10
#generate sparse matrix from count tables of patient P10 (initial and relapse)
total.5 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v005_total.TranscriptCounts.tsv.gz"))
total.6 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v006_total.TranscriptCounts.tsv.gz"))
total.17 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v017_total.TranscriptCounts.tsv.gz"))
total.18 <- as.matrix(read.table("2022-10-11_1st_batch/count_tables/EMB-JL-v018_total.TranscriptCounts.tsv.gz"))


total.5 <- as(total.5, "dgCMatrix")
total.6 <- as(total.6, "dgCMatrix")
total.17 <- as(total.17, "dgCMatrix")
total.18 <- as(total.18, "dgCMatrix")

total.73 <- as(total.73, "dgCMatrix")
total.74 <- as(total.74, "dgCMatrix")
total.77 <- as(total.77, "dgCMatrix")
total.78 <- as(total.78, "dgCMatrix")

ercc.total.5 <- grep("ERCC",rownames(total.5))
ercc.total.6 <- grep("ERCC",rownames(total.6))
ercc.total.17 <- grep("ERCC",rownames(total.17))
ercc.total.18 <- grep("ERCC",rownames(total.18))

#generate seurat objects for P10
total.object.5 <- CreateSeuratObject(counts = total.5[-ercc.total.5,], project = "EMB-JL-v005")
total.object.6 <- CreateSeuratObject(counts = total.6[-ercc.total.6,], project = "EMB-JL-v006")
total.object.17 <- CreateSeuratObject(counts = total.17[-ercc.total.17,], project = "EMB-JL-v017")
total.object.18 <- CreateSeuratObject(counts = total.18[-ercc.total.18,], project = "EMB-JL-v018")
#merge seurat objects
total.object.P10 <- merge(x = total.object.5, y = c(total.object.6, total.object.17, total.object.18), add.cell.ids = c("EMB-JL-v005", "EMB-JL-v006", "EMB-JL-v017", "EMB-JL-v018"))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P10-Ini', 384*2,1),matrix('P10-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P10)
total.object.P10 <- AddMetaData(
  object = total.object.P10,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)


#percentage human/mouse reads
total.object.P10[["human"]] <- PercentageFeatureSet(total.object.P10, pattern = "ENSG")
total.object.P10[["mouse"]] <- PercentageFeatureSet(total.object.P10, pattern = "ENSMUSG")
head(total.object.P10@meta.data, 5)
VlnPlot(total.object.P10, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P10 <- subset(total.object.P10, subset = human > 80 )
total.P10 <- subset(total.P10, features = rownames(total.P10)[grepl(rownames(total.P10),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total.P10$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total.P10$RNA) <- unique_row_names
head(total.P10@meta.data, 6)

#Quality control and selection of cells for further analysis
total.P10[["percent.mt"]] <- PercentageFeatureSet(total.object.P10, pattern = "MT")
head(total.P10@meta.data, 6)
VlnPlot(total.P10, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P10, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P10 <- subset(total.P10, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P10 <- NormalizeData(total.P10)
total.P10 <- FindVariableFeatures(total.P10, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P10), 10)
plot1 <- VariableFeaturePlot(total.P10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P10 <- ScaleData(total.P10)
total.P10 <- RunPCA(total.P10, features = VariableFeatures(object = total.P10))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P10 <- JoinLayers(total.P10)
total.P10 <- CellCycleScoring(total.P10, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P10 <- AddModuleScore(total.P10, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P10, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P10 <- ScaleData(total.P10, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P10))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P10.cellcycle <- RunPCA(total.P10, features = c(s.genes, g2m.genes))
DimPlot(total.P10.cellcycle, group.by = "Phase")
FeaturePlot(total.P10.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P10@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P10 <- RunPCA(object = total.P10, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P10), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P10, ndims = 50)
total.P10 <- FindNeighbors(total.P10, dims = 1:30)
total.P10 <- FindClusters(total.P10, resolution = 0.5)
total.P10 <- RunUMAP(total.P10, dims = 1:30)

#P1

#generate sparse matrix from count tables of patient P1 (initial and relapse)
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

total.103 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v103-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v103-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v103-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.104 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v104-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v104-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v104-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.107 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v107-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v107-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v107-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.108 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v108-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v108-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v108-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.103 <- grep("ERCC",rownames(total.103))
ercc.total.104 <- grep("ERCC",rownames(total.104))
ercc.total.107 <- grep("ERCC",rownames(total.107))
ercc.total.108 <- grep("ERCC",rownames(total.108))

#generate seurat objects for P1
total.object.103 <- CreateSeuratObject(counts = total.103[-ercc.total.103,], project = "EMB-JL-v103")
total.object.104 <- CreateSeuratObject(counts = total.104[-ercc.total.104,], project = "EMB-JL-v104")
total.object.107 <- CreateSeuratObject(counts = total.107[-ercc.total.107,], project = "EMB-JL-v107")
total.object.108 <- CreateSeuratObject(counts = total.108[-ercc.total.108,], project = "EMB-JL-v108")
#merge seurat objects
total.object.P1 <- merge(x = total.object.103, y = c(total.object.104, total.object.107, total.object.108), add.cell.ids = c("EMB-JL-v103", "EMB-JL-v104", "EMB-JL-v107", "EMB-JL-v108" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P1-Ini', 384*2,1),matrix('P1-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P1)
total.object.P1 <- AddMetaData(
  object = total.object.P1,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)


#percentage human/mouse reads
total.object.P1[["human"]] <- PercentageFeatureSet(total.object.P1, pattern = "ENSG")
total.object.P1[["mouse"]] <- PercentageFeatureSet(total.object.P1, pattern = "ENSMUSG")
head(total.object.P1@meta.data, 5)
VlnPlot(total.object.P1, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P1 <- subset(total.object.P1, subset = human > 80 )
total.P1 <- subset(total.P1, features = rownames(total.P1)[grepl(rownames(total.P1),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P1),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P1)<-gene_names

#Quality control and selection of cells for further analysis
total.P1[["percent.mt"]] <- PercentageFeatureSet(total.object.P1, pattern = "MT")
head(total.P1@meta.data, 6)
VlnPlot(total.P1, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P1 <- subset(total.P1, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P1 <- NormalizeData(total.P1)
total.P1 <- FindVariableFeatures(total.P1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P1), 10)
plot1 <- VariableFeaturePlot(total.P1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P1 <- ScaleData(total.P1)
total.P1 <- RunPCA(total.P1, features = VariableFeatures(object = total.P1))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P1 <- JoinLayers(total.P1)
total.P1 <- CellCycleScoring(total.P1, 
                              g2m.features = g2m.genes, 
                              s.features = s.genes)
total.P1 <- AddModuleScore(total.P1, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P1, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P1 <- ScaleData(total.P1, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P1))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P1.cellcycle <- RunPCA(total.P1, features = c(s.genes, g2m.genes))
DimPlot(total.P1.cellcycle, group.by = "Phase")
FeaturePlot(total.P1.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P1@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P1 <- RunPCA(object = total.P1, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P1), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P1, ndims = 50)
total.P1 <- FindNeighbors(total.P1, dims = 1:30)
total.P1 <- FindClusters(total.P1, resolution = 0.5)
total.P1 <- RunUMAP(total.P1, dims = 1:30)


#P3
#generate sparse matrix from count tables of patient P3 (initial and relapse)
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")
total.113 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v113-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v113-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v113-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.114 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v114-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v114-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v114-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.93 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v093-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v093-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v093-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.94 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v094-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v094-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v094-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)

ercc.total.113 <- grep("ERCC",rownames(total.113))
ercc.total.114 <- grep("ERCC",rownames(total.114))
ercc.total.93 <- grep("ERCC",rownames(total.93))
ercc.total.94 <- grep("ERCC",rownames(total.94))

#generate seurat objects for P3
total.object.113 <- CreateSeuratObject(counts = total.113[-ercc.total.113,], project = "EMB-JL-v113")
total.object.114 <- CreateSeuratObject(counts = total.114[-ercc.total.114,], project = "EMB-JL-v114")
total.object.93 <- CreateSeuratObject(counts = total.93[-ercc.total.93,], project = "EMB-JL-v093")
total.object.94 <- CreateSeuratObject(counts = total.94[-ercc.total.94,], project = "EMB-JL-v094")
#merge seurat objects
total.object.P3 <- merge(x = total.object.113, y = c(total.object.114, total.object.93, total.object.94), add.cell.ids = c("EMB-JL-v113", "EMB-JL-v114", "EMB-JL-v093", "EMB-JL-v094" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P3-Ini', 384*2,1),matrix('P3-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P3)
total.object.P3 <- AddMetaData(
  object = total.object.P3,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P3[["human"]] <- PercentageFeatureSet(total.object.P3, pattern = "ENSG")
total.object.P3[["mouse"]] <- PercentageFeatureSet(total.object.P3, pattern = "ENSMUSG")
head(total.object.P3@meta.data, 5)
VlnPlot(total.object.P3, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P3 <- subset(total.object.P3, subset = human > 80 )
total.P3 <- subset(total.P3, features = rownames(total.P3)[grepl(rownames(total.P3),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P3),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P3)<-gene_names

#Quality control and selection of cells for further analysis
total.P3[["percent.mt"]] <- PercentageFeatureSet(total.object.P3, pattern = "MT")
head(total.P3@meta.data, 6)
VlnPlot(total.P3, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P3 <- subset(total.P3, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P3 <- NormalizeData(total.P3)
total.P3 <- FindVariableFeatures(total.P3, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P3), 10)
plot1 <- VariableFeaturePlot(total.P3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P3 <- ScaleData(total.P3)
total.P3 <- RunPCA(total.P3, features = VariableFeatures(object = total.P3))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P3 <- JoinLayers(total.P3)
total.P3 <- CellCycleScoring(total.P3, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P3 <- AddModuleScore(total.P3, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P1, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P3 <- ScaleData(total.P3, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P3))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P3.cellcycle <- RunPCA(total.P3, features = c(s.genes, g2m.genes))
DimPlot(total.P3.cellcycle, group.by = "Phase")
FeaturePlot(total.P3.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P3@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P3 <- RunPCA(object = total.P3, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P3), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P3, ndims = 50)
total.P3 <- FindNeighbors(total.P3, dims = 1:30)
total.P3 <- FindClusters(total.P3, resolution = 0.5)
total.P3 <- RunUMAP(total.P3, dims = 1:30)

#P2
total.P2 <- readRDS("Analysis/P2_no_integration/total.P2.object.rds") #see Figure 1

#P12
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
#generate sparse matrix from count tables of patient P12 (initial and relapse)
total.9 <- as.matrix(read.table("9samples_human_mouse_mixed/raw/count_tables/EMB-JL-v009_total.TranscriptCounts.tsv.gz"))
total.10 <- as.matrix(read.table("9samples_human_mouse_mixed/raw/count_tables/EMB-JL-v010_total.TranscriptCounts.tsv.gz"))
total.21 <- as.matrix(read.table("9samples_human_mouse_mixed/raw/count_tables/EMB-JL-v021_total.TranscriptCounts.tsv.gz"))
total.22 <- as.matrix(read.table("9samples_human_mouse_mixed/raw/count_tables/EMB-JL-v022_total.TranscriptCounts.tsv.gz"))
total.9 <- as(total.9, "dgCMatrix")
total.10 <- as(total.10, "dgCMatrix")
total.21 <- as(total.21, "dgCMatrix")
total.22 <- as(total.22, "dgCMatrix")
ercc.total.9 <- grep("ERCC",rownames(total.9))
ercc.total.10 <- grep("ERCC",rownames(total.10))
ercc.total.21 <- grep("ERCC",rownames(total.21))
ercc.total.22 <- grep("ERCC",rownames(total.22))

#generate seurat objects
total.object.9 <- CreateSeuratObject(counts = total.9[-ercc.total.9,], project = "EMB-JL-v009")
total.object.10 <- CreateSeuratObject(counts = total.10[-ercc.total.10,], project = "EMB-JL-v010")
total.object.21 <- CreateSeuratObject(counts = total.21[-ercc.total.21,], project = "EMB-JL-v021")
total.object.22 <- CreateSeuratObject(counts = total.22[-ercc.total.22,], project = "EMB-JL-v022")
#merge seurat objects
total.object.P12 <- merge(x = total.object.9, y = c(total.object.10, total.object.21, total.object.22), add.cell.ids = c("EMB-JL-v009", "EMB-JL-v010", "EMB-JL-v021", "EMB-JL-v022" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P12-Ini', 384*2,1),matrix('P12-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P12)
total.object.P12 <- AddMetaData(
  object = total.object.P12,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P12[["human"]] <- PercentageFeatureSet(total.object.P12, pattern = "ENSG")
total.object.P12[["mouse"]] <- PercentageFeatureSet(total.object.P12, pattern = "ENSMUSG")
head(total.object.P12@meta.data, 5)
VlnPlot(total.object.P12, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P12 <- subset(total.object.P12, subset = human > 80 )
total.P12 <- subset(total.P12, features = rownames(total.P12)[grepl(rownames(total.P12),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total.P12$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total.P12$RNA) <- unique_row_names
head(total.P12@meta.data, 6)

#Quality control and selection of cells for further analysis
total.P12[["percent.mt"]] <- PercentageFeatureSet(total.object.P12, pattern = "MT")
head(total.P12@meta.data, 6)
VlnPlot(total.P12, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P12, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P12 <- subset(total.P12, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P12 <- NormalizeData(total.P12)
total.P12 <- FindVariableFeatures(total.P12, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P12), 10)
plot1 <- VariableFeaturePlot(total.P12)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P12 <- ScaleData(total.P12)
total.P12 <- RunPCA(total.P12, features = VariableFeatures(object = total.P12))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P12 <- JoinLayers(total.P12)
total.P12 <- CellCycleScoring(total.P12, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P12 <- AddModuleScore(total.P12, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P12, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P12 <- ScaleData(total.P12, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P12))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P12.cellcycle <- RunPCA(total.P12, features = c(s.genes, g2m.genes))
DimPlot(total.P2.cellcycle, group.by = "Phase")
FeaturePlot(total.P12.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P12@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P12 <- RunPCA(object = total.P12, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P12), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P12, ndims = 50)
total.P12 <- FindNeighbors(total.P12, dims = 1:30)
total.P12 <- FindClusters(total.P12, resolution = 0.5)
total.P12 <- RunUMAP(total.P12, dims = 1:30)


#P7
#generate sparse matrix from count tables of patient P7 (initial and relapse)
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")
total.85 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v085-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v085-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v085-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.86 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v086-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v086-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v086-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.97 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v097-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v097-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v097-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.98 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v098-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v098-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v098-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.85 <- grep("ERCC",rownames(total.85))
ercc.total.86 <- grep("ERCC",rownames(total.86))
ercc.total.97 <- grep("ERCC",rownames(total.97))
ercc.total.98 <- grep("ERCC",rownames(total.98))

#generate seurat objects for P7
total.object.85 <- CreateSeuratObject(counts = total.85[-ercc.total.85,], project = "EMB-JL-v085")
total.object.86 <- CreateSeuratObject(counts = total.86[-ercc.total.86,], project = "EMB-JL-v086")
total.object.97 <- CreateSeuratObject(counts = total.97[-ercc.total.97,], project = "EMB-JL-v097")
total.object.98 <- CreateSeuratObject(counts = total.98[-ercc.total.98,], project = "EMB-JL-v098")
#merge seurat objects
total.object.P7 <- merge(x = total.object.85, y = c(total.object.86, total.object.97, total.object.98), add.cell.ids = c("EMB-JL-v085", "EMB-JL-v086", "EMB-JL-v097", "EMB-JL-v098" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P7-Ini', 384*2,1),matrix('P7-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P7)
total.object.P7 <- AddMetaData(
  object = total.object.P7,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P7[["human"]] <- PercentageFeatureSet(total.object.P7, pattern = "ENSG")
total.object.P7[["mouse"]] <- PercentageFeatureSet(total.object.P7, pattern = "ENSMUSG")
head(total.object.P7@meta.data, 5)
VlnPlot(total.object.P7, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P7 <- subset(total.object.P7, subset = human > 80 )
total.P7 <- subset(total.P7, features = rownames(total.P7)[grepl(rownames(total.P7),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P7),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P7)<-gene_names

#Quality control and selection of cells for further analysis
total.P7[["percent.mt"]] <- PercentageFeatureSet(total.object.P7, pattern = "MT")
head(total.P7@meta.data, 6)
VlnPlot(total.P7, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P7, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P7 <- subset(total.P7, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P7 <- NormalizeData(total.P7)
total.P7 <- FindVariableFeatures(total.P7, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P7), 10)
plot1 <- VariableFeaturePlot(total.P7)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P7 <- ScaleData(total.P7)
total.P7 <- RunPCA(total.P7, features = VariableFeatures(object = total.P7))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P7 <- JoinLayers(total.P7)
total.P7 <- CellCycleScoring(total.P7, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P7 <- AddModuleScore(total.P7, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P7, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P7 <- ScaleData(total.P7, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P7))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P7.cellcycle <- RunPCA(total.P7, features = c(s.genes, g2m.genes))
DimPlot(total.P7.cellcycle, group.by = "Phase")
FeaturePlot(total.P7.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P7@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P7 <- RunPCA(object = total.P7, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P3), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P7, ndims = 50)
total.P7 <- FindNeighbors(total.P7, dims = 1:30)
total.P7 <- FindClusters(total.P7, resolution = 0.5)
total.P7 <- RunUMAP(total.P7, dims = 1:30)



##Fig 5a-c
#load object containing al patients with stemness marker values
obj.subgroups <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total_allsubgroups_after_integration_44k_genes.rds") ##see figure3 for scenic object containing all subgroups
#P6
# Extract stemness marker values for cells belonging to P6
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p6_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P6"]
stemness_values <- obj.subgroups@meta.data[p6_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P6 contains the same cells
matching_cells <- intersect(rownames(total.P6@meta.data), rownames(stemness_values))
# Transfer the values
total.P6@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P6_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P6, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#F0E442", "3" = "#56B4E9")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P6.png',width=3500,height=2500,res=600)
DimPlot(total.P6, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#F0E442", "3" = "#56B4E9")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P6.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P6, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#F0E442","#56B4E9" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()

#P8
# Extract stemness marker values for cells belonging to P8
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p8_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P8"]
stemness_values <- obj.subgroups@meta.data[p8_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P8 contains the same cells
matching_cells <- intersect(rownames(total.P8@meta.data), rownames(stemness_values))
# Transfer the values
total.P8@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P8_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P8, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P8.png',width=3500,height=2500,res=600)
DimPlot(total.P8, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P8.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P8, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                
dev.off()

#P10
# Extract stemness marker values for cells belonging to P10
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p10_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P10"]
stemness_values <- obj.subgroups@meta.data[p10_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P10 contains the same cells
matching_cells <- intersect(rownames(total.P10@meta.data), rownames(stemness_values))
# Transfer the values
total.P10@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P10_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P10, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P10.png',width=3500,height=2500,res=600)
DimPlot(total.P10, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P10.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P10, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P1
# Extract stemness marker values for cells belonging to P1
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p1_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P1"]
stemness_values <- obj.subgroups@meta.data[p1_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P1 contains the same cells
matching_cells <- intersect(rownames(total.P1@meta.data), rownames(stemness_values))
# Transfer the values
total.P1@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P1_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P1, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442", "4" = "#CC79A7", "5" = "#D55E00")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P1.png',width=3500,height=2500,res=600)
DimPlot(total.P1, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442", "4" = "#CC79A7", "5" = "#D55E00")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P1.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P1, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442", "#CC79A7","#D55E00" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P3
# Extract stemness marker values for cells belonging to P3
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p3_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P3"]
stemness_values <- obj.subgroups@meta.data[p3_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P3 contains the same cells
matching_cells <- intersect(rownames(total.P3@meta.data), rownames(stemness_values))
# Transfer the values
total.P3@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P3_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P3, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P3.png',width=3500,height=2500,res=600)
DimPlot(total.P3, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P3.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P3, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P2
# Extract stemness marker values for cells belonging to P2
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p2_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P2"]
stemness_values <- obj.subgroups@meta.data[p2_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P2 contains the same cells
matching_cells <- intersect(rownames(total.P2@meta.data), rownames(stemness_values))
# Transfer the values
total.P2@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P2_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P2, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9")) + NoLegend() + ggtitle("Stemness Score")  + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P2.png',width=3500,height=2500,res=600)
DimPlot(total.P2, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P2.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P2, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()   
dev.off()
#P12
# Extract stemness marker values for cells belonging to P12
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p12_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P12"]
stemness_values <- obj.subgroups@meta.data[p12_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P12 contains the same cells
matching_cells <- intersect(rownames(total.P12@meta.data), rownames(stemness_values))
# Transfer the values
total.P12@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P12_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P12, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P12.png',width=3500,height=2500,res=600)
DimPlot(total.P12, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P12.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P12, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend() 
dev.off()
#P7
# Extract stemness marker values for cells belonging to P7
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
p7_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P7"]
stemness_values <- obj.subgroups@meta.data[p7_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P12 contains the same cells
matching_cells <- intersect(rownames(total.P7@meta.data), rownames(stemness_values))
# Transfer the values
total.P7@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../Manuscript/figures/images/Fig3/Stemness_P7_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P7, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#56B4E9", "2" =  "#009E73", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2", "7" = "#666666", "8" = "#AD7700")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/UMAP_P7.png',width=3500,height=2500,res=600)
DimPlot(total.P7, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#56B4E9", "2" =  "#009E73", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2", "7" = "#666666", "8" = "#AD7700")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../Manuscript/figures/images/Fig3/barplot_P7.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P7, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00","#56B4E9", "#009E73","#F0E442", "#D55E00","#CC79A7", "#0072B2", "#666666", "#AD7700" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend() 
dev.off()

##Fig 5d
#calculate stem-like cluster enrichment of individual patients
Idents(total.P6) <- "seurat_clusters"
P6.stemlike <- subset(total.P6, subset = seurat_clusters == "3")
Idents(P6.stemlike) <- "Patient"
table(Idents(P6.stemlike))
#Ini: 0, Rel: 66
table(total.P6$Patient)#Ini: 586, Rel:618
0
66/618*100 = 10.68

Idents(total.P8) <- "seurat_clusters"
P8.stemlike <- subset(total.P8, subset = seurat_clusters == "2")
Idents(P8.stemlike) <- "Patient"
table(Idents(P8.stemlike))
#Ini: 3, Rel: 185
table(total.P8$Patient)#Ini: 646, Rel:630
3/646*100 = 0.46
185/630*100 = 29.37

Idents(total.P10) <- "seurat_clusters"
P10.stemlike <- subset(total.P10, subset = seurat_clusters == "2")
Idents(P10.stemlike) <- "Patient"
table(Idents(P10.stemlike))
#Ini: 5, Rel: 164
table(total.P10$Patient)#Ini: 626, Rel:704
5/626*100 = 0.80
164/704*100 = 23.3

Idents(total.P1) <- "seurat_clusters"
P1.stemlike <- subset(total.P1, subset = seurat_clusters == "2")
Idents(P1.stemlike) <- "Patient"
table(Idents(P1.stemlike))
#Ini: 9, Rel: 184
table(total.P1$Patient)#Ini: 416, Rel:626
9 / 416*100 = 2.16
184 / 626*100 = 29.39

Idents(total.P3) <- "seurat_clusters"
P3.stemlike <- subset(total.P3, subset = seurat_clusters == "2")
Idents(P3.stemlike) <- "Patient"
table(Idents(P3.stemlike))
#Ini: 4, Rel: 196
table(total.P3$Patient)#Ini: 702, Rel:653
4/702*100 = 0.57
196/653*100 = 30.02

Idents(total.P2) <- "seurat_clusters"
P2.stemlike <- subset(total.P2, subset = seurat_clusters == "2")
Idents(P2.stemlike) <- "Patient"
table(Idents(P2.stemlike))
#Ini: 9, Rel: 153
table(total.P2$Patient)#Ini: 653, Rel:578
9/653*100 = 1.38
153/578*100 = 26.47

Idents(total.P12) <- "seurat_clusters"
P12.stemlike <- subset(total.P12, subset = seurat_clusters == "2")
Idents(P12.stemlike) <- "Patient"
table(Idents(P12.stemlike))
#Ini: 0, Rel: 282
table(total.P12$Patient)#Ini: 650, Rel:648
0
282/648*100 = 43.52

Idents(total.P7) <- "seurat_clusters"
P7.stemlike <- subset(total.P7, subset = seurat_clusters == "1")
Idents(P7.stemlike) <- "Patient"
table(Idents(P7.stemlike))
#Ini: 12, Rel: 198
table(total.P7$Patient)#Ini: 700, Rel:669
12/700*100 = 1.71
198/669*100 = 29.60


Initial <- c(0, 0.46, 0.80, 2.16, 0.57, 1.38, 0, 1.71)
Relapse <- c(10.68, 29.37, 23.3, 29.39, 30.02, 26.47, 43.52, 29.6 )

my_data <- data.frame( 
  Disease_Stage = rep(c("Initial", "Relapse"), each = 8),
  Stem_like_Cluster = c( Initial, Relapse)
)

comparing_groups <- list(c("Initial", "Relapse") ) 
p <- ggpaired(my_data, x = "Disease_Stage", y = "Stem_like_Cluster", add = "dotplot", color = "Disease_Stage", palette =c("#56B4E9", "#56B4E9"), ylab = "Cluster Size in %", xlab = "Disease Stage")  
png('../../../../../../Manuscript/figures/images/Fig3/stem_like_cluster_barplot.png',width=4000,height=2500,res=600)
p + stat_compare_means(comparisons = comparing_groups, label.y = 35, method = "t.test", paired = TRUE, label = "p.format") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) + NoLegend() + ggtitle("Stem-like Cluster in Relapsing Patients")
dev.off()
#pval: 6.6e-05



