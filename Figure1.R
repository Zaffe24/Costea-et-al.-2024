###Figure 1: Expansion of a stem cell-like immature dormant subclone during the development of relapse in a child with T-ALL.###
library(readr)
library(stringr)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/")


##preprocessing Patient P2

#generate sparse matrix from count tables of patient P2 (initial and relapse)
total.45 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v045_total.TranscriptCounts.tsv.gz"))
total.46 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v046_total.TranscriptCounts.tsv.gz"))
total.53 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v053_total.TranscriptCounts.tsv.gz"))
total.54 <- as.matrix(read.table("2023-02-07_2nd_batch/count_tables/EMB-JL-v054_total.TranscriptCounts.tsv.gz"))
total.45 <- as(total.45, "dgCMatrix")
total.46 <- as(total.46, "dgCMatrix")
total.53 <- as(total.53, "dgCMatrix")
total.54 <- as(total.54, "dgCMatrix")
ercc.total.45 <- grep("ERCC",rownames(total.45))
ercc.total.46 <- grep("ERCC",rownames(total.46))
ercc.total.53 <- grep("ERCC",rownames(total.53))
ercc.total.54 <- grep("ERCC",rownames(total.54))

#generate seurat objects
total.object.45 <- CreateSeuratObject(counts = total.45[-ercc.total.45,], project = "EMB-JL-v045")
total.object.46 <- CreateSeuratObject(counts = total.46[-ercc.total.46,], project = "EMB-JL-v046")
total.object.53 <- CreateSeuratObject(counts = total.53[-ercc.total.53,], project = "EMB-JL-v053")
total.object.54 <- CreateSeuratObject(counts = total.54[-ercc.total.54,], project = "EMB-JL-v054")
#merge seurat objects
total.object.P2 <- merge(x = total.object.45, y = c(total.object.46, total.object.53, total.object.54), add.cell.ids = c("EMB-JL-v045", "EMB-JL-v046", "EMB-JL-v053", "EMB-JL-v054" ))

#add metadata with information on disease stage
cluster_letters.total <- as.character((rbind(matrix('P2-Ini', 384*2,1),matrix('P2-Rel', 384*2,1))))
names(cluster_letters.total) <- colnames(x=total.object.P2)
total.object.P2 <- AddMetaData(
  object = total.object.P2,
  metadata = cluster_letters.total,
  col.name = 'Patient'
)

#percentage human/mouse reads
total.object.P2[["human"]] <- PercentageFeatureSet(total.object.P2, pattern = "ENSG")
total.object.P2[["mouse"]] <- PercentageFeatureSet(total.object.P2, pattern = "ENSMUSG")
head(total.object.P2@meta.data, 5)
VlnPlot(total.object.P2, features = c("human", "mouse"), ncol = 2)
#remove mouse cells
total.P2 <- subset(total.object.P2, subset = human > 80 )
total.P2 <- subset(total.P2, features = rownames(total.P2)[grepl(rownames(total.P2),pattern = "ENSG")])

#keep only gene symbol
gene_names <- rownames(total.P2$RNA)
gene_names_modified <- sub(".*?-(.*?)-.*", "\\1", gene_names)
unique_row_names <- make.unique(gene_names_modified)
rownames(total.P2$RNA) <- unique_row_names
head(total.P2@meta.data, 6)

#Quality control and selection of cells for further analysis
total.P2[["percent.mt"]] <- PercentageFeatureSet(total.object.P2, pattern = "MT")
head(total.P2@meta.data, 6)
VlnPlot(total.P2, features = c("percent.mt"), ncol = 1)
VlnPlot(total.P2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
total.P2 <- subset(total.P2, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 5 & nCount_RNA < 100000)

#normalizing, scaling and PCA
total.P2 <- NormalizeData(total.P2)
total.P2 <- FindVariableFeatures(total.P2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(total.P2), 10)
plot1 <- VariableFeaturePlot(total.P2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
total.P2 <- ScaleData(total.P2)
total.P2 <- RunPCA(total.P2, features = VariableFeatures(object = total.P2))

#assign cell cycle scores based on S-phase genes, G2M-genes and canonical histone genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycling_histone_genes <- read_file("../canonical_histones.txt")
cycling_histone_genes <- read_lines(cycling_histone_genes)
cycling_histone_genes <- word(cycling_histone_genes, 2, sep="_")
total.P2 <- JoinLayers(total.P2)
total.P2 <- CellCycleScoring(total.P2, 
                             g2m.features = g2m.genes, 
                             s.features = s.genes)
total.P2 <- AddModuleScore(total.P2, features = list(cycling_histone_genes), name = "cycling_histones")

#PCA based on cell cycle scores: cells separate by their cell cycle phase prior to cell cycle regression
total.cell.cycle <- RunPCA(total.P2, features = c(s.genes, g2m.genes, cycling_histone_genes))
DimPlot(total.cell.cycle, group.by = "Phase")
FeaturePlot(total.cell.cycle, features = "cycling_histones1")

#regress out cell cycle scores during data scaling
total.P2 <- ScaleData(total.P2, vars.to.regress = c("S.Score", "G2M.Score", "cycling_histones1"), features = rownames(total.P2))

#When running a PCA based on cell cycle scores, cells no longer separate by cell cycle phase
total.P2.cellcycle <- RunPCA(total.P2, features = c(s.genes, g2m.genes))
DimPlot(total.P2.cellcycle, group.by = "Phase")
FeaturePlot(total.P2.cellcycle, features = "cycling_histones1")

#run PCA on regressed data and exclude noise genes + sex chromosomes 
all_genes <- rownames(total.P2@assays$RNA)
noise_genes <- all_genes[c(grep("^RPS", all_genes), grep("^RPL", all_genes),
                           grep("^AC0", all_genes), grep("^AL0", all_genes), grep("^AC1", all_genes),
                           grep("^AP0", all_genes), grep("^AL([0-9]+)", all_genes),
                           grep("^BX([0-9]+)", all_genes), grep("^C([0-9]+)", all_genes), grep("orf", all_genes))]
biomart <- readRDS("../msc_hg38_MAY2019_ageing_bioMart_2019_withSexChrom_full.RDS")
xy_genes <- unique(biomart[which(biomart$chromosome_name %in% c("chrY", "chrX")), "hgnc_symbol"])
exclude_genes <- unique(c(noise_genes, xy_genes))
total.P2 <- RunPCA(object = total.P2, verbose = FALSE, features = setdiff(VariableFeatures(object = total.P2), exclude_genes))

#cluster cells and run UMAP
ElbowPlot(total.P2, ndims = 50)
total.P2 <- FindNeighbors(total.P2, dims = 1:30)
total.P2 <- FindClusters(total.P2, resolution = 0.5)
total.P2 <- RunUMAP(total.P2, dims = 1:30)

setwd("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/2024-03-29_3rd_batch/2024-2184-emb-jl-v083-v114/")

##Fig 1b:UMAP initial and relapse
png('../../../../../../Manuscript/figures/images/Fig1/UMAP_P2.png',width=3500,height=2500,res=600)
DimPlot(total.P2, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9")) + ggtitle ("Cluster") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Fig 1d: UMAP initial, Fig 1e: UMAP relapse
P2.rel <-subset(total.P2, subset = Patient == 'P2-Rel')
P2.ini <-subset(total.P2, subset = Patient == 'P2-Ini')
png('../../../../../../Manuscript/figures/images/Fig1/UMAP_P2ini.png',width=3500,height=2500,res=600)
DimPlot(P2.ini, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9")) + ggtitle ("Initial Diagnosis") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()
png('../../../../../../Manuscript/figures/images/Fig1/UMAP_P2rel.png',width=3500,height=2500,res=600)
DimPlot(P2.rel, reduction = "umap", cols = c("0" = "#E69F00", "1" = "#009E73", "2" = "#56B4E9")) + ggtitle ("Relapse") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Fig 1c: barplot cluster frequency at initial and relapse
png('../../../../../../Manuscript/figures/images/Fig1/barplot_P2.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P2, "seurat_clusters", group.by = "Patient", color.panel = c("#E69F00", "#009E73","#56B4E9" )) + ylab("Cluster Frequency")+ xlab("Disease Stage") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold")) +  theme(text =element_text(size = 15)) + NoLegend()   
dev.off()

##Fig 1f: UMAP of predicted cell cycle phase based on cell cycle score
png('../../../../../../Manuscript/figures/images/Fig1/cell_cycle_umap_P2.png',width=3500,height=2500,res=600)
DimPlot(total.P2, reduction = "umap", group.by = "Phase", cols = c("G1" =  "#D99BBD", "G2M" = "#D5C711", "S" = "#A3A3A3"))+ ggtitle("Cell Cycle") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) 
dev.off()

##Fig 1g: Barplot of predicted cell cycle phase based on cell cycle score
png('../../../../../../Manuscript/figures/images/Fig1/cell_cycle_barplot_P2.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P2, "Phase",group.by = "seurat_clusters", color.panel = c("#D99BBD","#D5C711","#A3A3A3") ) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + NoLegend() + xlab("Seurat Clusters")
dev.off()

##Fig 1h: Cell type classification based on thymic single cell reference atlas (Jong-Eun Park et al, Science 2020)
#preprocessing thymic single cell reference
raw.data <- Read10X(data.dir = "../Annotation Single Cell Atlas T cell development Science/matrix_files/")
metadata <- read.csv("../Annotation Single Cell Atlas T cell development Science/metadata.csv")
head(metadata)
rownames(metadata) <- metadata$index
thymus <- CreateSeuratObject(counts = raw.data, meta.data = metadata, project = "ThymusReference")
thymus <- subset(x = thymus, downsample = 50000) 
thymus <- NormalizeData(thymus)
thymus <- FindVariableFeatures(thymus, selection.method = "vst", nfeatures = 2000)
thymus <- ScaleData(thymus)
thymus <- RunPCA(thymus, features = VariableFeatures(object = thymus))
ElbowPlot(thymus, ndims = 50)
thymus <- FindNeighbors(thymus, dims = 1:30)
thymus <- FindClusters(thymus, resolution = 0.5)
thymus <- RunUMAP(thymus, dims = 1:30, return.model = TRUE)
DimPlot(thymus, group.by = "Anno_level_fig1", label = TRUE)
saveRDS(thymus, file = "/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/Annotation Single Cell Atlas T cell development Science/thymus.rds")
#UMAP of cell type classification in patient P2
thymus <- readRDS("/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/Annotation Single Cell Atlas T cell development Science/thymus.rds")
anchor <- FindTransferAnchors(reference = thymus, query = total.P2,
                              dims = 1:30, reference.reduction = "pca")
total.P2 <- MapQuery(anchorset = anchor, reference = thymus, query = total.P2,
                     refdata = list(celltype = "Anno_level_fig1"), reference.reduction = "pca", reduction.model = "umap")
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_umap_P2.png',width=3500,height=2500,res=600)
DimPlot(total.P2, group.by = "predicted.celltype", reduction = "umap", cols = c("DN" = "#D55E00", "DP" = "#E69F00")) + ggtitle("Predicted Celltype") + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15))  
dev.off()

##Fig 1i: barplot of cell type classification in patient P2
png('../../../../../../Manuscript/figures/images/Fig1/cell_type_barplot_P2.png',width=3500,height=2500,res=600)
dittoBarPlot(total.P2, "predicted.celltype", group.by = "seurat_clusters", color.panel = c( "#D55E00",  "#E69F00" ))  + theme(plot.title = element_text(hjust = 0.5, face="bold"))  + theme(text =element_text(size = 15)) + ggtitle("") + NoLegend() + xlab("Seurat Clusters")
dev.off()

saveRDS(total.P2, "/g/korbel/Costea/Computational/VASAseq_data/TAL1_TALL/All_samples/Analysis/P2_no_integration/total.P2.object.rds")

##Fig 1j: regulons P2: output of pySCENIC pipeline 
#Add pySCENIC output to object and save it
scenic_df_wide <- read.csv("SCENIC/P2/new_aucell_mtx.tsv", 
                           sep = "\t", 
                           row.names = "Cell")
colnames(scenic_df_wide) <- colnames(scenic_df_wide) %>% str_replace(pattern = fixed("..."), "")
colnames(scenic_df_wide) <- colnames(scenic_df_wide) %>% str_replace(pattern = fixed("."), "-")
all_TFs <- colnames(scenic_df_wide)

total.P2 <- subset(total.P2,cells = rownames(scenic_df_wide))
names.use <- rownames(scenic_df_wide)[rownames(scenic_df_wide) %in% colnames(obj)]
scenic_df_wide <- scenic_df_wide[names.use, ]
total.P2[["scenic"]] <- CreateAssayObject(counts = t(scenic_df_wide))
#SCENIC output-target genes are assigned to each TF. Note: This is what we use for network re-construction
reg_files <- list.files("/g/scb2/zaugg/amathiou/2023July_TALL_Julia/SCENIC/P2/regulon", 
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
targene_df_raw <- reduce(df_list, bind_rows, .init = empty_df)
# make another copy with only the target genes that were used for the activity calculation 
targene_df <- filter(targene_df_raw, count > 100) #100 out of 250

#Marker genes
DefaultAssay(object = total.P2) <- "RNA"
DE_genes <- list()
for (p in unique(total.P2@meta.data$seurat_clusters)){
  tmp_de <- FindMarkers(object = total.P2, ident.1 = p)
  DE_genes[[paste0(p, "_+")]] <- rownames(tmp_de[which(tmp_de$avg_log2FC > 0.25 & tmp_de$p_val_adj < 0.000005), ])
  DE_genes[[paste0(p, "_-")]] <- rownames(tmp_de[which(tmp_de$avg_log2FC < -0.25 & tmp_de$p_val_adj < 0.000005), ])
}

pop_de <- list()
pop_de[[0]] <- unique(unlist(DE_genes[c("0_+", "0_-")]))
pop_de[[1]] <- unique(unlist(DE_genes[c("1_+", "1_-")]))
pop_de[[2]] <- unique(unlist(DE_genes[c("2_+", "2_-")]))

intersect(pop_de[[0]], unique(targene_df$target))
sig_TFs <- list()
fisher_df1 <- list()
fisher_df2 <- list()
for (p in c(0,1, 2)){
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
 
  #Plot enrichment
  
  q3 <- ggplot(table_tfs, aes(x = factor(TF, levels = levels_plot), y = log2OR, fill = padj)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          legend.position = "right") + xlab("TFs") + 
    ggtitle("Regulons of Cluster 2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(text = element_text(size = 20)) + 
    scale_fill_manual("padj", values = c("n.s." = "#BF812D", "<0.05" = "#35978F")) 
  
  png('../../../../../../Manuscript/figures/images/Fig1/regulons_barplot_P2.png',width=4500,height=2500,res=600)
  plot(q3)
  #dev.off()
  sig_TFs[[paste0(p, "_+")]] <- table_tfs[which(table_tfs$padj < 0.05 & table_tfs$log2OR > 0), "TF"]
  sig_TFs[[paste0(p, "_-")]] <- table_tfs[which(table_tfs$padj < 0.05 & table_tfs$log2OR < 0), "TF"]
  dev.off()
}


##Fig 1k: GSEA: see python script for analysis


