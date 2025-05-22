###Supplementary Figure 3

##preprocessing old libraries (TAL1 samples)
library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(clusterProfiler)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")

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

#total.73 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v073_total.TranscriptCounts.tsv.gz"))
#total.74 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v074_total.TranscriptCounts.tsv.gz"))
#total.77 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v077_total.TranscriptCounts.tsv.gz"))
#total.78 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v078_total.TranscriptCounts.tsv.gz"))

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


total.object <- merge(x = total.object.45, y = c(total.object.46, total.object.57, total.object.58, total.object.1, total.object.2, total.object.5, total.object.6, total.object.9, total.object.10, total.object.53, total.object.54, total.object.49, total.object.50, total.object.13, total.object.14, total.object.17, total.object.18, total.object.21, total.object.22, total.object.37, total.object.38, total.object.61, total.object.62, total.object.41, total.object.42, total.object.65, total.object.66, total.object.69, total.object.70), add.cell.ids = c("EMB-JL-v045", "EMB-JL-v046","EMB-JL-v057", "EMB-JL-v058", "EMB-JL-v001", "EMB-JL-v002", "EMB-JL-v005", "EMB-JL-v006", "EMB-JL-v009", "EMB-JL-v010", "EMB-JL-v053", "EMB-JL-v054", "EMB-JL-v049", "EMB-JL-v050", "EMB-JL-v013", "EMB-JL-v014", "EMB-JL-v017", "EMB-JL-v018", "EMB-JL-v021", "EMB-JL-v022", "EMB-JL-v037", "EMB-JL-v038","EMB-JL-v061", "EMB-JL-v062", "EMB-JL-v041", "EMB-JL-v042", "EMB-JL-v065", "EMB-JL-v066", "EMB-JL-v069", "EMB-JL-v070" ))

#add metadata with information on disease stage and patient of origin
cluster_letters.total <- as.character((rbind(matrix('P2-Ini', 384*2,1), matrix('P6-Ini', 384*2,1), matrix('P8-Ini', 384*2,1),matrix('P10-Ini', 384*2,1),matrix('P12-Ini', 384*2,1),matrix('P2-Rel', 384*2,1),matrix('P6-Rel', 384*2,1), matrix('P8-Rel', 384*2,1), matrix('P10-Rel', 384*2,1),matrix('P12-Rel', 384*2,1),matrix('P41', 384*2,1),matrix('P52', 384*2,1), matrix('P59', 384*2,1),matrix('P63', 384*2,1),matrix('P68', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#add metadata with information on disease group
cluster_letters.total <- as.character((rbind(matrix('Ini-R', 384*2,1), matrix('Ini-R', 384*2,1), matrix('Ini-R', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Rel', 384*2,1), matrix('Rel', 384*2,1), matrix('Rel', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-NR', 384*2,1),matrix('Ini-NR', 384*2,1), matrix('Ini-NR', 384*2,1),matrix('Ini-NR', 384*2,1),matrix('Ini-NR', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Disease_Group'
)

#add metadata with information on relapse type
cluster_letters.total <- as.character((rbind(matrix('type-2', 384*2,1), matrix('type-2', 384*2,1), matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1), matrix('type-2', 384*2,1), matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix(NA, 384*2,1),matrix(NA, 384*2,1), matrix(NA, 384*2,1),matrix(NA, 384*2,1),matrix(NA, 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Relapse_type'
)

#add metadata with information on T-ALL subgroup
cluster_letters.blast <- as.character((rbind(matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1))))
names(cluster_letters.blast) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.blast,
  col.name = 'Subgroup'
)

#add metadata with information on patient of origin
total.object@meta.data$Patients <- "placeholder"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P41")] <- "P41"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P52")] <- "P52"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P59")] <- "P59"
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

total.old.libraries <- total

##preprocessing new libraries (different T-ALL subgroups, type-1 and type-2 relapses)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

#P11
total.89 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v089-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v089-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v089-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.90 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v090-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v090-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v090-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.111 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v111-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v111-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v111-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.112 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v112-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v112-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v112-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)

ercc.total.89 <- grep("ERCC",rownames(total.89))
ercc.total.90 <- grep("ERCC",rownames(total.90))
ercc.total.111 <- grep("ERCC",rownames(total.111))
ercc.total.112 <- grep("ERCC",rownames(total.112))

total.object.89 <- CreateSeuratObject(counts = total.89[-ercc.total.89,], project = "EMB-JL-v089")
total.object.90 <- CreateSeuratObject(counts = total.90[-ercc.total.90,], project = "EMB-JL-v090")
total.object.111 <- CreateSeuratObject(counts = total.111[-ercc.total.111,], project = "EMB-JL-v111")
total.object.112 <- CreateSeuratObject(counts = total.112[-ercc.total.112,], project = "EMB-JL-v112")

#P1
total.103 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v103-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v103-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v103-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.104 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v104-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v104-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v104-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.107 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v107-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v107-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v107-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.108 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v108-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v108-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v108-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.103 <- grep("ERCC",rownames(total.103))
ercc.total.104 <- grep("ERCC",rownames(total.104))
ercc.total.107 <- grep("ERCC",rownames(total.107))
ercc.total.108 <- grep("ERCC",rownames(total.108))

total.object.103 <- CreateSeuratObject(counts = total.103[-ercc.total.103,], project = "EMB-JL-v103")
total.object.104 <- CreateSeuratObject(counts = total.104[-ercc.total.104,], project = "EMB-JL-v104")
total.object.107 <- CreateSeuratObject(counts = total.107[-ercc.total.107,], project = "EMB-JL-v107")
total.object.108 <- CreateSeuratObject(counts = total.108[-ercc.total.108,], project = "EMB-JL-v108")

#P5
total.83 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v083-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v083-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v083-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.84 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v084-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v084-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v084-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.95 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v095-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v095-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v095-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.96 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v096-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v096-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v096-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.83 <- grep("ERCC",rownames(total.83))
ercc.total.84 <- grep("ERCC",rownames(total.84))
ercc.total.95 <- grep("ERCC",rownames(total.95))
ercc.total.96 <- grep("ERCC",rownames(total.96))

total.object.83 <- CreateSeuratObject(counts = total.83[-ercc.total.83,], project = "EMB-JL-v083")
total.object.84 <- CreateSeuratObject(counts = total.84[-ercc.total.84,], project = "EMB-JL-v084")
total.object.95 <- CreateSeuratObject(counts = total.95[-ercc.total.95,], project = "EMB-JL-v095")
total.object.96 <- CreateSeuratObject(counts = total.96[-ercc.total.96,], project = "EMB-JL-v096")

#P3
total.113 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v113-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v113-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v113-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.114 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v114-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v114-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v114-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.93 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v093-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v093-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v093-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.94 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v094-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v094-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v094-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)

ercc.total.113 <- grep("ERCC",rownames(total.113))
ercc.total.114 <- grep("ERCC",rownames(total.114))
ercc.total.93 <- grep("ERCC",rownames(total.93))
ercc.total.94 <- grep("ERCC",rownames(total.94))

total.object.113 <- CreateSeuratObject(counts = total.113[-ercc.total.113,], project = "EMB-JL-v113")
total.object.114 <- CreateSeuratObject(counts = total.114[-ercc.total.114,], project = "EMB-JL-v114")
total.object.93 <- CreateSeuratObject(counts = total.93[-ercc.total.93,], project = "EMB-JL-v093")
total.object.94 <- CreateSeuratObject(counts = total.94[-ercc.total.94,], project = "EMB-JL-v094")

#P27
total.91 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v091-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v091-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v091-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.92 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v092-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v092-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v092-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.101 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v101-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v101-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v101-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.102 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v102-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v102-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v102-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.91 <- grep("ERCC",rownames(total.91))
ercc.total.92 <- grep("ERCC",rownames(total.92))
ercc.total.101 <- grep("ERCC",rownames(total.101))
ercc.total.102 <- grep("ERCC",rownames(total.102))

total.object.91 <- CreateSeuratObject(counts = total.91[-ercc.total.91,], project = "EMB-JL-v091")
total.object.92 <- CreateSeuratObject(counts = total.92[-ercc.total.92,], project = "EMB-JL-v092")
total.object.101 <- CreateSeuratObject(counts = total.101[-ercc.total.101,], project = "EMB-JL-v101")
total.object.102 <- CreateSeuratObject(counts = total.102[-ercc.total.102,], project = "EMB-JL-v102")

#P7
total.85 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v085-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v085-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v085-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.86 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v086-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v086-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v086-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.97 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v097-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v097-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v097-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.98 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v098-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v098-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v098-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.85 <- grep("ERCC",rownames(total.85))
ercc.total.86 <- grep("ERCC",rownames(total.86))
ercc.total.97 <- grep("ERCC",rownames(total.97))
ercc.total.98 <- grep("ERCC",rownames(total.98))

total.object.85 <- CreateSeuratObject(counts = total.85[-ercc.total.85,], project = "EMB-JL-v085")
total.object.86 <- CreateSeuratObject(counts = total.86[-ercc.total.86,], project = "EMB-JL-v086")
total.object.97 <- CreateSeuratObject(counts = total.97[-ercc.total.97,], project = "EMB-JL-v097")
total.object.98 <- CreateSeuratObject(counts = total.98[-ercc.total.98,], project = "EMB-JL-v098")

#P9
total.87 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v087-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v087-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v087-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.88 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v088-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v088-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v088-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.99 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v099-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v099-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v099-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.100 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v100-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v100-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v100-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.87 <- grep("ERCC",rownames(total.87))
ercc.total.88 <- grep("ERCC",rownames(total.88))
ercc.total.99 <- grep("ERCC",rownames(total.99))
ercc.total.100 <- grep("ERCC",rownames(total.100))

total.object.87 <- CreateSeuratObject(counts = total.87[-ercc.total.87,], project = "EMB-JL-v087")
total.object.88 <- CreateSeuratObject(counts = total.88[-ercc.total.88,], project = "EMB-JL-v088")
total.object.99 <- CreateSeuratObject(counts = total.99[-ercc.total.99,], project = "EMB-JL-v099")
total.object.100 <- CreateSeuratObject(counts = total.100[-ercc.total.100,], project = "EMB-JL-v100")

#P4
total.105 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v105-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v105-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v105-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.106 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v106-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v106-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v106-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.109 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v109-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v109-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v109-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
total.110 <- ReadMtx("raw_count_tables/poisson_corrected/EMB-JL-v110-raw/matrix.mtx", "raw_count_tables/poisson_corrected/EMB-JL-v110-raw/barcodes.tsv", "raw_count_tables/poisson_corrected/EMB-JL-v110-raw/features.tsv", cell.column = 1, feature.column = 1, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)


ercc.total.105 <- grep("ERCC",rownames(total.105))
ercc.total.106 <- grep("ERCC",rownames(total.106))
ercc.total.109 <- grep("ERCC",rownames(total.109))
ercc.total.110 <- grep("ERCC",rownames(total.110))

total.object.105 <- CreateSeuratObject(counts = total.105[-ercc.total.105,], project = "EMB-JL-v105")
total.object.106 <- CreateSeuratObject(counts = total.106[-ercc.total.106,], project = "EMB-JL-v106")
total.object.109 <- CreateSeuratObject(counts = total.109[-ercc.total.109,], project = "EMB-JL-v109")
total.object.110 <- CreateSeuratObject(counts = total.110[-ercc.total.110,], project = "EMB-JL-v110")

total.object <- merge(x = total.object.89, y = c( total.object.90, total.object.111, total.object.112, total.object.103, total.object.104, total.object.107, total.object.108, total.object.83, total.object.84, total.object.95, total.object.96, total.object.113, total.object.114, total.object.93, total.object.94, total.object.91, total.object.92, total.object.101, total.object.102, total.object.85, total.object.86, total.object.97, total.object.98, total.object.87, total.object.88, total.object.99, total.object.100, total.object.105, total.object.106, total.object.109, total.object.110), add.cell.ids = c( "EMB-JL-v089", "EMB-JL-v090", "EMB-JL-v111", "EMB-JL-v112", "EMB-JL-v103", "EMB-JL-v104", "EMB-JL-v107", "EMB-JL-v108", "EMB-JL-v083", "EMB-JL-v084", "EMB-JL-v095", "EMB-JLv096", "EMB-JL-v113", "EMB-JL-v114", "EMB-JL-v093", "EMB-JL-v094", "EMB-JL-v091", "EMB-JL-v092", "EMB-JL-v101", "EMB-JL-v102", "EMB-JL-v085", "EMB-JL-v086", "EMB-JL-v097", "EMB-JL-v098", "EMB-JL-v087", "EMB-JL-v088", "EMB-JL-v099", "EMB-JL-v100", "EMB-JL-v105", "EMB-JL-v106", "EMB-JL-v109", "EMB-JL-v110" ))

#add metadata with information on disease stage and patient of origin
cluster_letters.total <- as.character((rbind(matrix('P11-Ini', 384*2,1),matrix('P11-Rel', 384*2,1),matrix('P1-Ini', 384*2,1),matrix('P1-Rel', 384*2,1),matrix('P5-Ini', 384*2,1),matrix('P5-Rel', 384*2,1),matrix('P3-Ini', 384*2,1),matrix('P3-Rel', 384*2,1),matrix('P27-Ini', 384*2,1),matrix('P27-Rel', 384*2,1),matrix('P7-Ini', 384*2,1),matrix('P7-Rel', 384*2,1),matrix('P9-Ini', 384*2,1),matrix('P9-Rel', 384*2,1),matrix('P4-Ini', 384*2,1),matrix('P4-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#add metadata with information on disease group
cluster_letters.total <- as.character((rbind(matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1),matrix('Ini-R', 384*2,1),matrix('Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Disease_Group'
)

#add metatdata with information on relapse type
cluster_letters.total <- as.character((rbind(matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-2', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1),matrix('type-1', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Relapse_type'
)

#add metadata with information on T-ALL subgroup
cluster_letters.total <- as.character((rbind(matrix("TAL1", 384*2,1),matrix("TAL1", 384*2,1),matrix("NKX2", 384*2,1),matrix("NKX2", 384*2,1),matrix("NKX2", 384*2,1),matrix("NKX2", 384*2,1),matrix("HOXA", 384*2,1),matrix("HOXA", 384*2,1),matrix('HOXA', 384*2,1),matrix('HOXA', 384*2,1),matrix('TLX1', 384*2,1),matrix('TLX1', 384*2,1),matrix('TLX1', 384*2,1),matrix('TLX1', 384*2,1),matrix('TLX1', 384*2,1),matrix('TLX1', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object)
total.object <- AddMetaData(
  object = total.object,
  metadata = cluster_letters.total,
  col.name = 'Subgroup'
)

#add metadata with information on patient of origin
total.object@meta.data$Patients <- "placeholder"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P11-Ini")] <- "P11"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P11-Rel")] <- "P11"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P1-Ini")] <- "P1"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P1-Rel")] <- "P1"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P5-Ini")] <- "P5"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P5-Rel")] <- "P5"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P3-Ini")] <- "P3"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P3-Rel")] <- "P3"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P27-Ini")] <- "P27"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P27-Rel")] <- "P27"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P7-Ini")] <- "P7"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P7-Rel")] <- "P7"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P9-Ini")] <- "P9"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P9-Rel")] <- "P9"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P4-Ini")] <- "P4"
total.object@meta.data$Patients[which(total.object@meta.data$Patient == "P4-Rel")] <- "P4"

#Percentage human/mouse reads
total.object[["human"]] <- PercentageFeatureSet(total.object, pattern = "ENSG")
total.object[["mouse"]] <- PercentageFeatureSet(total.object, pattern = "ENSMUSG")
head(total.object@meta.data, 5)
VlnPlot(total.object, features = c("human", "mouse"), ncol = 2, group.by = "Patient")
organism.plot <- FeatureScatter(total.object, feature1 = "human", feature2 = "mouse")
organism.plot

#remove mouse cells
total <- subset(total.object, subset = human > 80 )
total <- subset(total, features = rownames(total)[grepl(rownames(total),pattern = "ENSG")])

#keep only gene symbols
gene_names <- rownames(total$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total$RNA) <- unique_row_names

total.new.libraries <- total

##merge at this step old and new libraries from preprocessing
total <- merge(x = total.old.libraries, y = total.new.libraries)
head(total)

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
total <- FindClusters(total, resolution = 0.5, cluster.name = "unintegrated_clusters")
total <- RunUMAP(total, dims = 1:50, reduction.name = "umap.unintegrated")

#CCA integration
total.integrated <- IntegrateLayers(object = total, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                    verbose = FALSE)
total.integrated[["RNA"]] <- JoinLayers(total.integrated[["RNA"]])
total.integrated <- FindNeighbors(total.integrated, reduction = "integrated.cca", dims = 1:50)
total.integrated <- FindClusters(total.integrated, resolution = 0.5)
total.integrated <- RunUMAP(total.integrated, dims = 1:50, reduction = "integrated.cca")

##Suppl Fig 3a: UMAP patients
png('../../../../../../Manuscript/figures/images/Fig3/patients_umapintegrated_allsubgroups.png',width=3500,height=2500,res=600)
DimPlot(total.integrated, group.by = "Patients", reduction = "umap") + ggtitle("All Patients after CCA Integration") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))
dev.off()

##Suppl Fig 3b: UMAP clusters
png('../../../../../../Manuscript/figures/images/Fig3/clusters_umapintegrated_allsubgroups.png',width=3500,height=2500,res=600)
DimPlot(total.integrated, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#56B4E9", "2" =  "#009E73", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2", "7" = "#666666", "8" = "#AD7700", "9" = "#FF61CC")) + ggtitle("Clusters after CCA Integration") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))
dev.off()

##Suppl Fig 3c: UMAP stemness score
png('../../../../../../Manuscript/figures/images/Fig3/stemness_umapintegrated_allsubgroups.png',width=3500,height=2500,res=600)
FeaturePlot(total.integrated, reduction = "umap", features = "stemness_marker1") + ggtitle("Stemness Score") + scale_color_viridis_c() + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Suppl Fig 3d: barplot
png('../../../../../../Manuscript/figures/images/Fig3/vlnplot_umapintegrated_stemness.png',width=3500,height=2500,res=600)
VlnPlot(total.integrated, features= "stemness_marker1", group.by="seurat_clusters", pt.size = 0, cols = c("0" = "#E69F00", "1" = "#56B4E9", "2" =  "#009E73", "3" = "#F0E442", "4" = "#D55E00", "5" ="#CC79A7", "6" ="#0072B2", "7" = "#666666", "8" = "#AD7700", "9" = "#FF61CC")) + NoLegend() + ggtitle("Stemness Score") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +xlab("Cluster") + theme(text =element_text(size = 15)) & geom_hline(yintercept = 0.25, col = 'red', linetype = "dashed") 
dev.off()

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

combined_plot <- wrap_plots(plots, ncol = 2) + 
  plot_annotation(
    title = "Enrichment of TAL1 top10 stem-like markers in cells of all subgroups",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

png('../../../../../Manuscript/figures/images/Fig3/DGEtop10_flipped.png',width=7000,height=11000,res=600)
print(combined_plot)
dev.off()


##Suppl Fig 3d:
DefaultAssay(obj.subgroups) <- "scenic"
Idents(obj.subgroups) <- "Patients"
features = c("FOSL2", "ETS1", "KLF2", "NFKB1", "ZMIZ1", "FLI1", "KLF6", "FOS", "FOXO1", "RELB", "JUN", "NFKB2", "FOSB", "E2F3", "E2F2", "BRCA1", "E2F1", "TFDP1")
png('../../../../../../Manuscript/figures/images/Fig3/regulon_per_patient_dotplot.png',width=6500,height=5500,res=600)
DotPlot(obj.subgroups, features = features, dot.scale = 10) + scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00") + RotatedAxis() + ggtitle("Regulon Activity per Patient") + labs(colour="Average Regulon Activity") +  theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 20)) & scale_y_discrete(limit=c("P41","P63", "P6","P8", "P10",  "P1", "P3", "P59", "P2", "P12", "P7", "P52", "P68", "P11",  "P5",  "P27", "P9", "P4")) &xlab("TFs")
dev.off()



