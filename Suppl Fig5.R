library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)



setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

##preprocessing individual patients with stem-like cells < 5%
#P11
#generate sparse matrix from count tables of patient P11 (initial and relapse)
total.89 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v089-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v089-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v089-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.90 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v090-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v090-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v090-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.111 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v111-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v111-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v111-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.112 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v112-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v112-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v112-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.89 <- grep("ERCC",rownames(total.89))
ercc.total.90 <- grep("ERCC",rownames(total.90))
ercc.total.111 <- grep("ERCC",rownames(total.111))
ercc.total.112 <- grep("ERCC",rownames(total.112))

#generate seurat object for P11
total.object.89 <- CreateSeuratObject(counts = total.89[-ercc.total.89,], project = "EMB-JL-v089")
total.object.90 <- CreateSeuratObject(counts = total.90[-ercc.total.90,], project = "EMB-JL-v090")
total.object.111 <- CreateSeuratObject(counts = total.111[-ercc.total.111,], project = "EMB-JL-v111")
total.object.112 <- CreateSeuratObject(counts = total.112[-ercc.total.112,], project = "EMB-JL-v112")
#merge seurat objects
total.object.P11 <- merge(x = total.object.89, y = c(total.object.90, total.object.111, total.object.112), add.cell.ids = c("EMB-JL-v089", "EMB-JL-v090", "EMB-JL-v111", "EMB-JL-v112" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P11-Ini', 384*2,1),matrix('P11-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P11)
total.object.P11 <- AddMetaData(
  object = total.object.P11,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P11[["human"]] <- PercentageFeatureSet(total.object.P11, pattern = "ENSG")
total.object.P11[["mouse"]] <- PercentageFeatureSet(total.object.P11, pattern = "ENSMUSG")
head(total.object.P11@meta.data, 5)
VlnPlot(total.object.P11, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P11 <- subset(total.object.P11, subset = human > 80 )
total.P11 <- subset(total.P11, features = rownames(total.P11)[grepl(rownames(total.P11),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P11),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P11)<-gene_names

#Quality control and selection of cells for further analysis
total.P11[["percent.mt"]] <- PercentageFeatureSet(total.object.P11, pattern = "MT")
head(total.P11@meta.data, 6)
VlnPlot(total.P11, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P11, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P11 <- subset(total.P11, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P11 <- NormalizeData(total.P11)
total.P11 <- FindVariableFeatures(total.P11, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P11), 10)
plot1 <- VariableFeaturePlot(total.P11)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P11 <- ScaleData(total.P11)
total.P11 <- RunPCA(total.P11, features = VariableFeatures(object = total.P11))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P11 <- JoinLayers(total.P11)
total.P11 <- CellCycleScoring(total.P11, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P11 <- AddModuleScore(total.P11, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P11, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P11 <- ScaleData(total.P11, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P11))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P11.cellcycle <- RunPCA(total.P11, features = c(s.genes, g2m.genes))
DimPlot(total.P11.cellcycle, group.by = "Phase")
FeaturePlot(total.P11.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P11@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P11 <- RunPCA(object = total.P11, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P11), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P11, ndims = 50)
total.P11 <- FindNeighbors(total.P11, dims = 1:30)
total.P11 <- FindClusters(total.P11, resolution = 0.5)
total.P11 <- RunUMAP(total.P11, dims = 1:30)

#P5
#generate sparse matrix from count tables of patient P5 (initial and relapse)
total.83 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v083-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v083-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.84 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v084-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v084-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v084-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.95 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v095-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v095-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v095-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.96 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v096-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v096-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v096-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.83 <- grep("ERCC",rownames(total.83))
ercc.total.84 <- grep("ERCC",rownames(total.84))
ercc.total.95 <- grep("ERCC",rownames(total.95))
ercc.total.96 <- grep("ERCC",rownames(total.96))
#generate seurat objects for P5
total.object.83 <- CreateSeuratObject(counts = total.83[-ercc.total.83,], project = "EMB-JL-v083")
total.object.84 <- CreateSeuratObject(counts = total.84[-ercc.total.84,], project = "EMB-JL-v084")
total.object.95 <- CreateSeuratObject(counts = total.95[-ercc.total.95,], project = "EMB-JL-v095")
total.object.96 <- CreateSeuratObject(counts = total.96[-ercc.total.96,], project = "EMB-JL-v096")
#merge seurat objects
total.object.P5 <- merge(x = total.object.83, y = c(total.object.84, total.object.95, total.object.96), add.cell.ids = c("EMB-JL-v083", "EMB-JL-v084", "EMB-JL-v095", "EMB-JL-v096" ))
#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P5-Ini', 384*2,1),matrix('P5-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P5)
total.object.P5 <- AddMetaData(
  object = total.object.P5,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P5[["human"]] <- PercentageFeatureSet(total.object.P5, pattern = "ENSG")
total.object.P5[["mouse"]] <- PercentageFeatureSet(total.object.P5, pattern = "ENSMUSG")
head(total.object.P5@meta.data, 5)
VlnPlot(total.object.P5, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P5 <- subset(total.object.P5, subset = human > 80 )
total.P5 <- subset(total.P5, features = rownames(total.P5)[grepl(rownames(total.P5),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P5),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P5)<-gene_names

#Quality control and selection of cells for further analysis
total.P5[["percent.mt"]] <- PercentageFeatureSet(total.object.P5, pattern = "MT")
head(total.P5@meta.data, 6)
VlnPlot(total.P5, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P5, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P5 <- subset(total.P5, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P5 <- NormalizeData(total.P5)
total.P5 <- FindVariableFeatures(total.P5, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P5), 10)
plot1 <- VariableFeaturePlot(total.P5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P5 <- ScaleData(total.P5)
total.P5 <- RunPCA(total.P5, features = VariableFeatures(object = total.P5))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P5 <- JoinLayers(total.P5)
total.P5 <- CellCycleScoring(total.P5, 
                              g2m.features = g2m.genes, 
                              s.features = s.genes)
total.P5 <- AddModuleScore(total.P5, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P5, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P5 <- ScaleData(total.P5, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P5))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P5.cellcycle <- RunPCA(total.P5, features = c(s.genes, g2m.genes))
DimPlot(total.P5.cellcycle, group.by = "Phase")
FeaturePlot(total.P5.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P5@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P5 <- RunPCA(object = total.P5, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P5), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P5, ndims = 50)
total.P5 <- FindNeighbors(total.P5, dims = 1:30)
total.P5 <- FindClusters(total.P5, resolution = 0.5)
total.P5 <- RunUMAP(total.P5, dims = 1:30)

#P27
#generate sparse matrix from count tables of patient P27 (initial and relapse)
total.91 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v091-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v091-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v091-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.92 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v092-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v092-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v092-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.101 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v101-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v101-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v101-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.102 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v102-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v102-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v102-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.91 <- grep("ERCC",rownames(total.91))
ercc.total.92 <- grep("ERCC",rownames(total.92))
ercc.total.101 <- grep("ERCC",rownames(total.101))
ercc.total.102 <- grep("ERCC",rownames(total.102))
#generate seurat objects
total.object.91 <- CreateSeuratObject(counts = total.91[-ercc.total.91,], project = "EMB-JL-v091")
total.object.92 <- CreateSeuratObject(counts = total.92[-ercc.total.92,], project = "EMB-JL-v092")
total.object.101 <- CreateSeuratObject(counts = total.101[-ercc.total.101,], project = "EMB-JL-v101")
total.object.102 <- CreateSeuratObject(counts = total.102[-ercc.total.102,], project = "EMB-JL-v102")
#merge seurat objects
total.object.P27 <- merge(x = total.object.91, y = c(total.object.92, total.object.101, total.object.102), add.cell.ids = c("EMB-JL-v091", "EMB-JL-v092", "EMB-JL-v101", "EMB-JL-v102" ))
#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P27-Ini', 384*2,1),matrix('P27-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P27)
total.object.P27 <- AddMetaData(
  object = total.object.P27,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P27[["human"]] <- PercentageFeatureSet(total.object.P27, pattern = "ENSG")
total.object.P27[["mouse"]] <- PercentageFeatureSet(total.object.P27, pattern = "ENSMUSG")
head(total.object.P27@meta.data, 5)
VlnPlot(total.object.P27, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P27 <- subset(total.object.P27, subset = human > 80 )
total.P27 <- subset(total.P27, features = rownames(total.P27)[grepl(rownames(total.P27),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P27),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P27)<-gene_names

#Quality control and selection of cells for further analysis
total.P27[["percent.mt"]] <- PercentageFeatureSet(total.object.P27, pattern = "MT")
head(total.P27@meta.data, 6)
VlnPlot(total.P27, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P27, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P27 <- subset(total.P27, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P27 <- NormalizeData(total.P27)
total.P27 <- FindVariableFeatures(total.P27, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P27), 10)
plot1 <- VariableFeaturePlot(total.P27)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P27 <- ScaleData(total.P27)
total.P27 <- RunPCA(total.P27, features = VariableFeatures(object = total.P27))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P27 <- JoinLayers(total.P27)
total.P27 <- CellCycleScoring(total.P27, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P27 <- AddModuleScore(total.P27, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P27, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P27 <- ScaleData(total.P27, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P27))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P27.cellcycle <- RunPCA(total.P27, features = c(s.genes, g2m.genes))
DimPlot(total.P27.cellcycle, group.by = "Phase")
FeaturePlot(total.P27.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P27@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P27 <- RunPCA(object = total.P27, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P27), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P27, ndims = 50)
total.P27 <- FindNeighbors(total.P27, dims = 1:30)
total.P27 <- FindClusters(total.P27, resolution = 0.5)
total.P27 <- RunUMAP(total.P27, dims = 1:30)

#P9
#generate sparse matrix from count tables of patient P9 (initial and relapse)
total.87 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v087-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v087-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v087-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.88 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v088-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v088-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v088-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.99 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v099-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v099-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v099-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.100 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v100-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v100-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v100-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.87 <- grep("ERCC",rownames(total.87))
ercc.total.88 <- grep("ERCC",rownames(total.88))
ercc.total.99 <- grep("ERCC",rownames(total.99))
ercc.total.100 <- grep("ERCC",rownames(total.100))
#generate seurat objects
total.object.87 <- CreateSeuratObject(counts = total.87[-ercc.total.87,], project = "EMB-JL-v087")
total.object.88 <- CreateSeuratObject(counts = total.88[-ercc.total.88,], project = "EMB-JL-v088")
total.object.99 <- CreateSeuratObject(counts = total.99[-ercc.total.99,], project = "EMB-JL-v099")
total.object.100 <- CreateSeuratObject(counts = total.100[-ercc.total.100,], project = "EMB-JL-v100")
#merge seurat objects
total.object.P9 <- merge(x = total.object.87, y = c(total.object.88, total.object.99, total.object.100), add.cell.ids = c("EMB-JL-v87", "EMB-JL-v088", "EMB-JL-v099", "EMB-JL-v100" ))
#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P9-Ini', 384*2,1),matrix('P9-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P9)
total.object.P9 <- AddMetaData(
  object = total.object.P9,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P9[["human"]] <- PercentageFeatureSet(total.object.P9, pattern = "ENSG")
total.object.P9[["mouse"]] <- PercentageFeatureSet(total.object.P9, pattern = "ENSMUSG")
head(total.object.P9@meta.data, 5)
VlnPlot(total.object.P9, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P9 <- subset(total.object.P9, subset = human > 80 )
total.P9 <- subset(total.P9, features = rownames(total.P9)[grepl(rownames(total.P9),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P9),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P9)<-gene_names

#Quality control and selection of cells for further analysis
total.P9[["percent.mt"]] <- PercentageFeatureSet(total.object.P9, pattern = "MT")
head(total.P9@meta.data, 6)
VlnPlot(total.P9, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P9, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P9 <- subset(total.P9, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P9 <- NormalizeData(total.P9)
total.P9 <- FindVariableFeatures(total.P9, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P9), 10)
plot1 <- VariableFeaturePlot(total.P9)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P9 <- ScaleData(total.P9)
total.P9 <- RunPCA(total.P9, features = VariableFeatures(object = total.P9))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P9 <- JoinLayers(total.P9)
total.P9 <- CellCycleScoring(total.P9, 
                              g2m.features = g2m.genes, 
                              s.features = s.genes)
total.P9 <- AddModuleScore(total.P9, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P9, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P9 <- ScaleData(total.P9, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P9))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P9.cellcycle <- RunPCA(total.P9, features = c(s.genes, g2m.genes))
DimPlot(total.P9.cellcycle, group.by = "Phase")
FeaturePlot(total.P9.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P9@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P9 <- RunPCA(object = total.P9, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P9), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P9, ndims = 50)
total.P9 <- FindNeighbors(total.P9, dims = 1:30)
total.P9 <- FindClusters(total.P9, resolution = 0.5)
total.P9 <- RunUMAP(total.P9, dims = 1:30)

#P4
#generate sparse matrix from count tables of patient P4 (initial and relapse)
total.105 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v105-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v105-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v105-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.106 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v106-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v106-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v106-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.109 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v109-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v109-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v109-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.110 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v110-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v110-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v110-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
ercc.total.105 <- grep("ERCC",rownames(total.105))
ercc.total.106 <- grep("ERCC",rownames(total.106))
ercc.total.109 <- grep("ERCC",rownames(total.109))
ercc.total.110 <- grep("ERCC",rownames(total.110))
#generate seurat objects
total.object.105 <- CreateSeuratObject(counts = total.105[-ercc.total.105,], project = "EMB-JL-v105")
total.object.106 <- CreateSeuratObject(counts = total.106[-ercc.total.106,], project = "EMB-JL-v106")
total.object.109 <- CreateSeuratObject(counts = total.109[-ercc.total.109,], project = "EMB-JL-v109")
total.object.110 <- CreateSeuratObject(counts = total.110[-ercc.total.110,], project = "EMB-JL-v110")
#merge seurat objects
total.object.P4 <- merge(x = total.object.105, y = c(total.object.106, total.object.109, total.object.110), add.cell.ids = c("EMB-JL-v105", "EMB-JL-v106", "EMB-JL-v109", "EMB-JL-v110" ))
#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P4-Ini', 384*2,1),matrix('P4-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P4)
total.object.P4 <- AddMetaData(
  object = total.object.P4,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P4[["human"]] <- PercentageFeatureSet(total.object.P4, pattern = "ENSG")
total.object.P4[["mouse"]] <- PercentageFeatureSet(total.object.P4, pattern = "ENSMUSG")
head(total.object.P4@meta.data, 5)
VlnPlot(total.object.P4, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P4 <- subset(total.object.P4, subset = human > 80 )
total.P4 <- subset(total.P4, features = rownames(total.P4)[grepl(rownames(total.P4),pattern = "ENSG")])

#change to gene names
table <- read_tsv("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", col_names = FALSE)
table<-data.frame(table)
rownames(table)<-table$X1
gene_names <- table[rownames(total.P4),'X2']
gene_names<-make.unique(gene_names,sep = '.')
rownames(total.P4)<-gene_names

#Quality control and selection of cells for further analysis
total.P4[["percent.mt"]] <- PercentageFeatureSet(total.object.P4, pattern = "MT")
head(total.P4@meta.data, 6)
VlnPlot(total.P4, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P4, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P4 <- subset(total.P4, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P4 <- NormalizeData(total.P4)
total.P4 <- FindVariableFeatures(total.P4, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P4), 10)
plot1 <- VariableFeaturePlot(total.P4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P4 <- ScaleData(total.P4)
total.P4 <- RunPCA(total.P4, features = VariableFeatures(object = total.P4))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P4 <- JoinLayers(total.P4)
total.P4 <- CellCycleScoring(total.P4, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P4 <- AddModuleScore(total.P4, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P4, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P4 <- ScaleData(total.P4, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P4))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P4.cellcycle <- RunPCA(total.P4, features = c(s.genes, g2m.genes))
DimPlot(total.P4.cellcycle, group.by = "Phase")
FeaturePlot(total.P4.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P4@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P4 <- RunPCA(object = total.P4, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P4), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P4, ndims = 50)
total.P4 <- FindNeighbors(total.P4, dims = 1:30)
total.P4 <- FindClusters(total.P4, resolution = 0.5)
total.P4 <- RunUMAP(total.P4, dims = 1:30)

##Suppl Fig 5a-c: umap, vlnplot and boxplot of patients with <5% stem like cells
#load object containing al patients with stemness marker values
obj.subgroups <- readRDS("/g/korbel/Costea/Computational/SCENIC/2023July_TALL_Julia/total_allsubgroups_after_integration_44k_genes.rds") ##see figure3 for scenic object containing all subgroups
#P11
# Extract stemness marker values for cells belonging to P11
p11_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P11"]
stemness_values <- obj.subgroups@meta.data[p11_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P27 contains the same cells
matching_cells <- intersect(rownames(total.P11@meta.data), rownames(stemness_values))
# Transfer the values
total.P11@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../../../Manuscript/figures/images/Fig4/Stemness_P11_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P11, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P11.png',width=3500,height=2500,res=600)
DimPlot(total.P11, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/barplot_P11.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P11, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P5
# Extract stemness marker values for cells belonging to P5
p5_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P5"]
stemness_values <- obj.subgroups@meta.data[p5_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P5 contains the same cells
matching_cells <- intersect(rownames(total.P5@meta.data), rownames(stemness_values))
# Transfer the values
total.P5@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../../../Manuscript/figures/images/Fig4/Stemness_P5_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P5, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#666666", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + ylim(-0.15, 0.3) + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P5.png',width=3500,height=2500,res=600)
DimPlot(total.P5, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#666666", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/barplot_P5.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P5, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73", "#666666",  "#F0E442",  "#D55E00", "#CC79A7", "#0072B2")) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P27
# Extract stemness marker values for cells belonging to P27
p27_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P27"]
stemness_values <- obj.subgroups@meta.data[p27_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P27 contains the same cells
matching_cells <- intersect(rownames(total.P27@meta.data), rownames(stemness_values))
# Transfer the values
total.P27@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../../../Manuscript/figures/images/Fig4/Stemness_P27_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P27, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#0077BE", "1" = "#009E73", "2" = "#E69F00" , "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7","6" = "#666666" )) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P27.png',width=3500,height=2500,res=600)
DimPlot(total.P27, reduction = "umap", cols = c("0" = "#0077BE", "1" = "#009E73", "2" = "#E69F00" , "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" = "#666666")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/barplot_P27.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P27, "seurat_clusters", group.by = "Patient", color.panel = c("#0077BE", "#009E73","#E69F00","#F0E442", "#D55E00","#CC79A7", "#666666")) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P9
# Extract stemness marker values for cells belonging to P9
p9_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P9"]
stemness_values <- obj.subgroups@meta.data[p9_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P27 contains the same cells
matching_cells <- intersect(rownames(total.P9@meta.data), rownames(stemness_values))
# Transfer the values
total.P9@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../../../Manuscript/figures/images/Fig4/Stemness_P9_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P9, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00" , "1" = "#009E73", "2" =  "#666666", "3" = "#F0E442","4" = "#D55E00", "5" ="#CC79A7")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P9.png',width=3500,height=2500,res=600)
DimPlot(total.P9, reduction = "umap", cols = c("0" = "#E69F00" , "1" = "#009E73", "2" =  "#666666", "3" = "#F0E442","4" = "#D55E00", "5" ="#CC79A7")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/barplot_P9.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P9, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00" ,  "#009E73",  "#666666", "#F0E442", "#D55E00", "#CC79A7")) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()
#P4
# Extract stemness marker values for cells belonging to P4
p4_cells <- rownames(obj.subgroups@meta.data)[obj.subgroups@meta.data$Patients == "P4"]
stemness_values <- obj.subgroups@meta.data[p4_cells, "stemness_marker1", drop = FALSE]
# Ensure total.P27 contains the same cells
matching_cells <- intersect(rownames(total.P4@meta.data), rownames(stemness_values))
# Transfer the values
total.P4@meta.data[matching_cells, "stemness_marker1"] <- stemness_values[matching_cells, , drop = FALSE]
#generate plots
png('../../../../../../Manuscript/figures/images/Fig4/Stemness_P4_vlnplot.png',width=3500,height=2500,res=600)
VlnPlot(total.P4, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#0077BE", "3" = "#F0E442")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P4.png',width=3500,height=2500,res=600)
DimPlot(total.P4, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#0077BE", "3" = "#F0E442")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig4/barplot_P4.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P4, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#0077BE","#F0E442" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()                                                                                                                                                                                                                                  
dev.off()

#Suppl Fig5d: visualization of technical replicates to exclude batch effects leading to cluster composition differences in initial vs relapse
#relapsing patients with > 5% stem-like cells (for P2 check Suppl Fig1)
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")
total.P6 <- readRDS("Analysis/P6_no_integration/total.P6.object_swapped.rds")
total.P8 <- readRDS("Analysis/P8_no_integration/total.P8.object.rds")
total.P10 <- readRDS(file = "Analysis/P10_no_integration/total.P10_PDXonly.rds")
total.P2 <- readRDS("Analysis/P2_no_integration/total.P2.object.rds")
total.P12 <- readRDS("Analysis/P12_no_integration/P12.object")
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")
total.P1 <- readRDS("../../Analysis/P1_no_integration/P1.object.rds")
total.P3 <- readRDS("../../Analysis/P3_no_integration/P3.object.rds")
total.P7 <- readRDS("../../Analysis/P7_no_integration/P7.object.rds")

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P6replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P6, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P6 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P8replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P8, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P8 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P10replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P10, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P10 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P1replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P1, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P1 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P3replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P3, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P3 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P12replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P12, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P12 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P7replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P7, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P7 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##technical replicates of relapsing patients with <5% stem like cells
setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")
total.P11 <- readRDS("../../Analysis/P11_no_integration/P11.object.rds")
total.P5 <- readRDS("../../Analysis/P5_no_integration/P5.object.rds")
total.P27 <- readRDS("../../Analysis/P27_no_integration/P27.object.rds")
total.P9 <- readRDS("../../Analysis/P9_no_integration/P9.object.rds")
total.P4 <- readRDS("../../Analysis/P4_no_integration/P4.object.rds")

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P27replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P27, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P27 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P9replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P9, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P9 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P4replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P4, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P4 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

#batch effect visible
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P11replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P11, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P11 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

#one plate with fewer cells
png('../../../../../../Manuscript/figures/images/Fig4/UMAP_P5replicates.png',width=4000,height=2500,res=600)
DimPlot(total.P5, reduction = "umap", group.by = "orig.ident") + ggtitle ("Replicates of P5 Initial and Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

