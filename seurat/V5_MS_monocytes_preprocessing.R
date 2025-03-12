# Calling libraries
# Seurat v5
library(openxlsx)
library(Seurat)
library(harmony)
library(reticulate)
library(BPCells)
library(dplyr)
library(SeuratWrappers)
library(Azimuth)
library(Matrix)
library(sctransform)
library(car)
library(scater)
library(ggplot2)
library(patchwork)
library(ggrepel)
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")
library(BiocParallel)
register(MulticoreParam(14))
library(SeuratDisk)
library(velocyto.R)
library(ggpubr)



#########################         BPCells           ###################

## Set working dir
setwd("/path/to/working/dir/")

################### First round of sequence
## Loop through h5 files and output BPCells matrices on-disk
file.dir <- "../../cellranger/"
files.set <- c(
  "GEM1",
  "GEM10",
  "GEM11",
  "GEM12",
  "GEM13",
  "GEM14",
  "GEM15",
  "GEM16",
  "GEM17",
  "GEM18",
  "GEM19",
  "GEM2",
  "GEM20",
  "GEM21",
  "GEM22",
  "GEM23",
  "GEM24",
  "GEM25",
  "GEM3",
  "GEM4",
  "GEM5",
  "GEM6",
  "GEM7",
  "GEM8",
  "GEM9"
)

data.list <- c()

for (i in 1:length(files.set)) {
  ## Create binary matrices
  path <- paste0(file.dir, files.set[i], "/outs/filtered_feature_bc_matrix.h5")
  data <- open_matrix_10x_hdf5(path = path)
  write_matrix_dir(mat = data, dir = paste0(files.set[i], "_BP2"))
  ## Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP2"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  dataset_name <- files.set[i]
  mat2 <- CreateSeuratObject(counts = mat, min.features = 500) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    subset(subset = nFeature_RNA < 5000 & percent.mt < 10
    ) %>% NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")
  mat2$GEM <- dataset_name
  data.list[[i]] <- mat2
  rm(mat)
  rm(mat2)
  rm(data)
}

#################### Second round of sequence

file.dir <- "../../new_MS_mono_objs/seurat_objs/"
files.set <- c(
  "GEM26",
  "GEM27",
  "GEM28",
  "GEM29",
  "GEM30"
)
nums <- c("1","2","3","4","5")

for (i in 1:length(files.set)) {
  ## Create binary matrices
  path <- paste0(file.dir, "hashtag", nums[i],".rds")
  data <- readRDS(path)
  data[["RNA"]]$counts <- as(object = data[["RNA"]]$counts, Class = "dgCMatrix")
  write_matrix_dir(mat = data[["RNA"]]$counts, dir = paste0(files.set[i], "_BP"))
  ## Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP"))
  #mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  dataset_name <- files.set[i]
  mat2 <- CreateSeuratObject(counts = mat, min.features = 500, min.cells = 3) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    subset(subset = nFeature_RNA < 5000 & percent.mt < 10 
    ) %>% NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")
  mat2$GEM <- dataset_name
  data.list[[i+25]] <- mat2
  rm(mat)
  rm(mat2)
  rm(data)
}

####### Merge Both data

files.set <- c(
  "GEM1",
  "GEM10",
  "GEM11",
  "GEM12",
  "GEM13",
  "GEM14",
  "GEM15",
  "GEM16",
  "GEM17",
  "GEM18",
  "GEM19",
  "GEM2",
  "GEM20",
  "GEM21",
  "GEM22",
  "GEM23",
  "GEM24",
  "GEM25",
  "GEM3",
  "GEM4",
  "GEM5",
  "GEM6",
  "GEM7",
  "GEM8",
  "GEM9",
  "GEM26",
  "GEM27",
  "GEM28",
  "GEM29",
  "GEM30"
)

# Name layers
names(data.list) <- files.set

# Merge layers and create seurat obj during merging
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 1500)
prot.combined <- merge(data.list[[1]], y = data.list[2:length(data.list)],
                       add.cell.ids = files.set, merge.data = T)
VariableFeatures(prot.combined) <- features

# Normalize and scale merged obj
prot.combined <- NormalizeData(prot.combined)
prot.combined <- FindVariableFeatures(prot.combined)
prot.combined <- ScaleData(prot.combined)

## Dimensionality reduction and integration
prot.combined <- RunPCA(prot.combined)
gc()
ElbowPlot(prot.combined, ndims = 50)
prot.combined <- FindNeighbors(prot.combined, dims = 1:40, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 0.5, cluster.name = "unintegrated_clusters")

prot.combined <- RunUMAP(prot.combined, dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

## add metadata
meta.data <- read.csv("meta.txt", header = T)
prot.combined$prog <- prot.combined$GEM
prot.combined$disease <- prot.combined$GEM
prot.combined$sex <- prot.combined$GEM
prot.combined$pair <- prot.combined$GEM
prot.combined$age <- prot.combined$GEM
prot.combined$edss <- prot.combined$GEM
prot.combined$disease.state <- prot.combined$GEM

for (i in 1:length(meta.data$GEM)) {
  prot.combined$prog <- recode(prot.combined$prog, "meta.data$GEM[i] = meta.data$Progression[i]")
  prot.combined$sex <- recode(prot.combined$sex, "meta.data$GEM[i] = meta.data$Sex[i]")
  prot.combined$disease <- recode(prot.combined$disease, "meta.data$GEM[i] = meta.data$Disease[i]")
  prot.combined$pair <- recode(prot.combined$pair, "meta.data$GEM[i] = meta.data$Pair[i]")
  prot.combined$age <- recode(prot.combined$age, "meta.data$GEM[i] = meta.data$Age[i]")
  prot.combined$edss <- recode(prot.combined$edss, "meta.data$GEM[i] = meta.data$EDSS[i]")
  prot.combined$disease.state <- recode(prot.combined$disease.state, "meta.data$GEM[i] = meta.data$Group[i]")
}
prot.combined$prog_sex <- paste(prot.combined$sex, prot.combined$prog, sep = "_")
prot.combined$disease_sex <- paste(prot.combined$sex, prot.combined$disease, sep = "_")
prot.combined$disease.state_sex <- paste(prot.combined$sex, prot.combined$disease.state, sep = "_")


## Save and Load data
saveRDS(object = prot.combined, file = "obj_BP2_unintegrated_annotated.Rds")
prot.combined <- readRDS("./obj_BP2_unintegrated_annotated.Rds")


##########

prot.combined <- JoinLayers(prot.combined)

myeloids <- subset(prot.combined, idents = c("8","2","6","0","1","3","4","5","9","17"))
myeloids <- RenameIdents(myeloids, `0` = "Classical", `1` = "Classical", `3` = "Classical",
                         `2` = "Nonclassical", `6` = "Intermediate", `4` = "Classical",
                         `5` = "Classical", `9` = "Classical",`17` = "Classical",`8` = "mo-DC")
myeloids$celltypes <- Idents(myeloids)

## Save and Load data
saveRDS(object = myeloids, file = "myeloid_BP2_unintegrated_annotated2.Rds")
## Reload
# myeloids <- readRDS("./myeloid_BP2_unintegrated_annotated2.Rds")


######### ############    celltypes records

library(ggstatsplot)
myeloids$disease.patient <- paste(myeloids$disease.state, myeloids$GEM, sep = "_")
num.cells <- as.data.frame(table(myeloids$disease.patient, myeloids$disease.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$disease.patient, myeloids$seurat_clusters))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.cluster_disease_patient.csv")


########## cell types frequencies

df_freqs <- read.csv("freq.num.cells.cluster_disease_patient.csv", header = T)
celltypes2 <- colnames(df_freqs)[11:35]
coluna1 <- c("HC","RRMS","Non-progressor","Progressor")

my_comparisons <- list(c("HC","RRMS"),
                       c("HC","Non-progressor"),
                       c("HC","Progressor"),
                       c("Non-progressor","RRMS"),
                       c("Progressor","RRMS"),
                       c("Progressor","Non-progressor"))

## do it by cluster
p1 <- ggplot(df_freqs, aes_string(x="factor(progression,levels=coluna1)", y="X4", color="progression")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + # Add pairwise comparisons p-value
  theme_minimal() + rotate_x_text(angle = 45) + ylab("Cell percentage (%)") + 
  theme(axis.text.x = element_text(face="bold", size=12),axis.text.y = element_text(face="bold", size=12)) +
  ggtitle("Cluster 4") +
  theme(plot.title = element_text(hjust = 0.5))
#theme_classic()
p1

###################################### extract embeddings and clusters

umap_coords <- Embeddings(myeloids, "umap.unintegrated")
rownames(umap_coords) <- gsub("_",":",rownames(umap_coords))
rownames(umap_coords) <- gsub("-1","x",rownames(umap_coords))
colnames(umap_coords) <- c("umap_1","umap_2")
umap_coords1 <- as.matrix(umap_coords[grep("GEM10:",rownames(umap_coords)),])

cluster_data <- data.frame(cell = colnames(myeloids), seurat_clusters = myeloids$seurat_clusters)
cluster_data$cell <- gsub("_",":",cluster_data$cell)
rownames(cluster_data) <- gsub("_",":",rownames(cluster_data))
cluster_data$cell <- gsub("-1","x",cluster_data$cell)
rownames(cluster_data) <- gsub("-1","x",rownames(cluster_data))
bm <- AddMetaData(object = bm,metadata = cluster_data, col.name = "seurat_clusters")
Idents(bm) <- "seurat_clusters"

cluster_data <- data.frame(cell = rownames(cluster_data), seurat_clusters = cluster_data$seurat_clusters)
cluster_data1 <- as.matrix(cluster_data[grep("GEM10:",cluster_data$seurat_clusters),])
cluster_data1 <- cluster_data[rownames(cluster_data) %in% rownames(cluster_data)[grep("GEM10:",rownames(cluster_data))],]


## Save V5 as V3
myeloids[["RNA"]] <- JoinLayers(myeloids[["RNA"]])
myeloids[["RNA"]]$scale.data <- NULL
myeloids[["RNA3"]] <- as(object = myeloids[["RNA"]], Class = "Assay")
myeloids2 <- myeloids

myeloids2[["RNA"]] <- myeloids2[["RNA3"]]
myeloids2[["RNA3"]] <- NULL

## Save and Load data
saveRDS(object = myeloids2, file = "myeloid_BP2_unintegrated_annotated2.v3.Rds")
# rm(myeloids)
# rm(myeloids2)
myeloids2 <- readRDS("./myeloid_BP2_unintegrated_annotated.v3.Rds")

## rename clusters
seurat_clusters <- c("0","1","2","3","4","5","6","8","9","17")
new_clusters <- c("0","1","2","3","4","5","6","7","8","9")
myeloids$new_clusters <- myeloids$seurat_clusters

for (i in 1:length(new_clusters)) {
  myeloids$new_clusters <- recode(myeloids$new_clusters, "seurat_clusters[i] = new_clusters[i]")
}


### Diff exp analysis

myeloids$celltypes.prog <- paste(myeloids$celltypes, myeloids$prog, sep = "_")
Idents(myeloids) <- "celltypes.prog"

cell.types <- c("Nonclassical","Intermediate","Classical","mo-DC")
diseases <- c("_HC","_non-progressor","_RR")

#myeloids <- JoinLayers(myeloids)
for (i in 1:length(cell.types)) {
  for (e in diseases) {
    zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_progressor")
    ,ident.2 = paste0(cell.types[i], e)
    , slot = "data",
    assay = "RNA",
    features = NULL,
    logfc.threshold = 0.0,
    test.use = "wilcox",
    min.pct = 0.4,
    min.diff.pct = -Inf,
    verbose = TRUE,
    only.pos = FALSE,
    max.cells.per.ident = Inf,
    random.seed = 1,
    latent.vars = NULL,
    min.cells.feature = 3,
    min.cells.group = 3,
    pseudocount.use = 1,
    mean.fxn = NULL,
    fc.name = NULL,
    base = 2,
    densify = FALSE,
    recorrect_umi = TRUE
    )
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    zk.response0 <- zk.response0[order(zk.response0$avg_log2FC),]
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_prog_progressor_x",e,"_", cell.types[i], "_DEGs_minpct0.4.xlsx"))
    rm(zk.response0) 
  }
}

myeloids$cluster.prog <- paste(myeloids$seurat_clusters, myeloids$prog, sep = "_")
Idents(myeloids) <- "cluster.prog"
nums <- c("0","1","2","3","4","5","6","8","9","17")

for (i in 1:length(nums)) {
  for (e in diseases) {
    zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(nums[i], "_progressor")
                                ,ident.2 = paste0(nums[i], e)
                                , slot = "data",
                                assay = "RNA",
                                features = NULL,
                                logfc.threshold = 0.0,
                                test.use = "wilcox",
                                min.pct = 0.4,
                                min.diff.pct = -Inf,
                                verbose = TRUE,
                                only.pos = FALSE,
                                max.cells.per.ident = Inf,
                                random.seed = 1,
                                latent.vars = NULL,
                                min.cells.feature = 3,
                                min.cells.group = 3,
                                pseudocount.use = 1,
                                mean.fxn = NULL,
                                fc.name = NULL,
                                base = 2,
                                densify = FALSE,
                                recorrect_umi = TRUE
    )
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    zk.response0 <- zk.response0[order(zk.response0$avg_log2FC),]
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_prog_progressor_x",e,"_", nums[i], "_DEGs_minpct0.4.xlsx"))
    rm(zk.response0) 
  }
}


