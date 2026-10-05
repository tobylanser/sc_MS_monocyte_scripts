# Seurat v5
library(openxlsx)
library(MAST)
library(Seurat)
library(harmony)
library(reticulate)
library(BPCells)
library(dplyr)
library(RColorBrewer)
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
library(paletteer)
library(nichenetr)
library(tidyverse)
library(circlize)
library(ggpubr)
library(VennDiagram)
library(pheatmap)
library(gplots)
library(GOplot)
library(readr)
library(robustbase)
library(EnhancedVolcano)
library(alluvial)
library(reshape2)
library(scico)
library(ggsci)
library(rcartocolor)
library(ggside)
library(viridis)
library(ggstatsplot)
library(grid)
library(shadowtext)
library(tidyr)



#########################         BPCells           ###################

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_jakel_el.al_2019/Seurat/")

expr_dt <- data.table::fread(file = "../preprocessed/MSCtr_snRNA_ExpressionMatrix_R.txt")
expr_mat <- as.data.frame(expr_dt)
rownames(expr_mat) <- expr_mat[[1]]
expr_mat[[1]] <- NULL
expr_mat <- as.matrix(expr_mat)
rm(expr_dt)
gc()

meta <- read.delim("../preprocessed/MSCtr_snRNA_FinalAnnotationTable.txt", row.names = 1, header = T)

prot.combined <- CreateSeuratObject(counts = expr_mat, min.cells = 1, min.features = 100, meta.data = meta) %>% 
  PercentageFeatureSet(pattern = "^MT.", col.name = "percent.mt") %>% 
  subset(subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 10 #& nCount_RNA < 40000
  ) %>% NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
# prot.combined <- AddMetaData(object = prot.combined, metadata = meta)

# meta.data <- read.delim("meta.txt", header = T)
rm(expr_mat)
gc()

# Visualize QC metrics as a violin plot
# prot.combined$cohort <- "jakel"
# VlnPlot(prot.combined, pt.size = 0.0,group.by = "cohort",raster = F,
#        features = c("nFeature_RNA", "nCount_RNA", "MT.ND4"), ncol = 3)
# plot1 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster = F, group.by = "cohort")
# plot2 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F, group.by = "cohort")
# plot1 + plot2

### Normalize and scale merged obj
# prot.combined <- NormalizeData(prot.combined, normalization.method = "LogNormalize")
# prot.combined <- FindVariableFeatures(prot.combined, selection.method = "vst", nfeatures = 3000)
prot.combined <- ScaleData(prot.combined, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                           model.use = "linear")#, features = all.genes)
gc()
### Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(prot.combined), 10)

### plot variable features with and without labels
#plot1 <- VariableFeaturePlot(prot.combined, raster = F)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, raster = F)
#plot1 + plot2

### Dimensionality reduction and integration
prot.combined <- RunPCA(prot.combined, npcs = 50)
gc()
ElbowPlot(prot.combined, ndims = 50)
prot.combined <- FindNeighbors(prot.combined, dims = 1:40, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 0.5, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        label = T,
        repel = T,
        #group.by = "seurat_clusters"
        # ncol = 2,
        # split.by = "pathology"
)#, combine = F)

## Save and Load V5 data
saveRDS(object = prot.combined, file = "obj_unintegrated_scaled_by.MT.counts.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_scaled_by.MT.counts.Rds")



## cell types annotation

DotPlot(prot.combined, features = c( "GFAP","SLC1A3","AQP4","LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4",   #Astrocytes
                                     "SLC17A6","SLC17A7","NRGN","CAMK2A", "SATB2", "COL5A1","SDK2","NEFM","HTR2C",    #Excitatory_neurons
                                     "SLC32A1","GAD1","GAD2","TAC1","PENK","SST","NPY","MYBPC1","PVALB","GABBR2",   #Inhibitory_neurons
                                     "OLIG2", "MBP","MOBP","PLP1","MOG","CLDN11","MYRF","GALC","ERMN","MAG",   #Oligodendrocytes
                                     "VCAN","CSPG4","PDGFRA", "SOX10","NEU4", "PCDH15","GPR37L1","C1QL1","CDO1","EPN2",   #Oligodendrocyte_precursor_cells
                                     "AMBP","HIGD1B","COX4I2", "AOC3","PDE5A","PTH1R","P2RY14","ABCC9","KCNJ8","CD248",  #Pericytes
                                     "FLT1","CLDN5", "VTN","ITM2A", "VWF", "FAM167B","BMX","CLEC1B",    #Endothelial_cells
                                     "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",   #Microglia
                                     "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2"#, #Macrophages
                                     # "CD8A","CD4","CD3E","CD3D","CD19","CD22","PTPRC","CD27"
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = T, 
# group.by = "type_broad",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()


prot.combined[["RNA"]] <- JoinLayers(prot.combined[["RNA"]])
prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap.unintegrated", reorder = T)
prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = T, min.pct = 0.40, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.clusters.markers.xlsx")
prot.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5


DotPlot(prot.combined, features = unique(top5$gene),
        col.max = 20, 
        dot.scale = 10, 
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

Idents(prot.combined) <- "unintegrated_clusters"

prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Oligodendrocytes",
                              `1`= "Astrocytes",                     # Homeostatic Microglia
                              `2`= "Excitatory neurons",  
                              `3`= "Astrocytes",
                              `4`= "Inhibitory neurons", 
                              `5`= "Oligodendrocytes",
                              `6`= "Excitatory neurons",
                              `7`= "Inhibitory neurons",
                              `8`= "Excitatory neurons",
                              `9`= "Oligodendrocytes",
                              `10`= "Oligodendrocyte precursor cells",
                              `11`= "Excitatory neurons",
                              `12`= "Oligodendrocytes",
                              `13`= "Excitatory neurons",
                              `14`= "Inhibitory neurons",                   # 
                              `15`= "Excitatory neurons",
                              `16`= "Microglia",
                              `17`= "Inhibitory neurons",
                              `18`= "Microglia",
                              `19`= "Excitatory neurons",
                              `20`= "Excitatory neurons",
                              `21`= "Excitatory neurons",
                              `22`= "Excitatory neurons",
                              `23`= "Oligodendrocyte precursor cells",
                              `24`= "Endothelial cells",
                              `25`= "Oligodendrocytes",
                              `26`= "Inhibitory neurons",          # 
                              `27`= "Astrocytes",              # 
                              `28`= "Excitatory neurons",
                              `29`= "Astrocytes", # Astrocytes
                              `30`= "Oligodendrocytes",           #
                              `31`= "Excitatory neurons",
                              `32`= "Excitatory neurons",           # 
                              `33`= "Oligodendrocytes",         
                              `34`= "Microglia", # Foam cells
                              `35`= "Inhibitory neurons",
                              `36`= "Oligodendrocytes",
                              `37`= "Oligodendrocyte precursor cells",
                              `38`= "Microglia",             
                              `39`= "Astrocytes",                      # 
                              `40`= "Excitatory neurons",
                              `41`= "Excitatory neurons",                  # 
                              `42`= "Excitatory neurons",
                              `43`= "Inhibitory neurons",
                              `44`= "Excitatory neurons",       
                              `45`= "Neuroendocrine cells" 
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"

## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.annotated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats.annotated.Rds")

##### freq
num.cells <- as.data.frame(table(prot.combined$pathology, prot.combined$pathology))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(prot.combined$pathology, prot.combined$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.prot.combined_disease_lesion.csv")

paletteLength <- 9
colors <- colorRampPalette( rev(brewer.pal(11, "Set3")))(paletteLength)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
barplot(tfreq.num.cells.celltype, col = colors,#sample(col_vector, n), 
        legend.text = rownames(tfreq.num.cells.celltype),
        xlim = c(0,7.2),
        main = "Cell types abundance", xlab = "", ylab = "Cells percentage (%)", cex.names = 0.9)

myeloids$lesion.patient <- paste(myeloids$lesion_type, myeloids$GEM, sep = "_")
num.cells <- as.data.frame(table(myeloids$lesion.patient, myeloids$lesion.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$lesion.patient, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.myeloids_lesion_patient.csv")




clus4.markers <- c("SDC2","PPIF","PTGS2","NLRP3","HBEGF","IL1B","NFKB1","FOSL2","GPR183","EGR1","EGR2","EGR3","CXCL8",
                   "G0S2","RGCC","TRIB1")

prot.combined <- AddModuleScore(
  prot.combined,
  features = clus4.markers,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "clus4.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)










gene <- "SDC2"
FeaturePlot(prot.combined, features = gene, 
            raster = F,
            # pt.size = 2.,
            # ncol = 3,
            reduction = "umap.unintegrated", 
            # min.cutoff = 0.7,
            # max.cutoff = 1,
            cols = c("gray90","#9E1021"),
            # split.by = "Lesion"
) + theme(legend.position = "right")

# Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "type_broad"
prot.combined$pathology <- factor(prot.combined$pathology, levels = c("control WM","control cortex","myelinated cortex","demyelinated cortex"))

my_comparisons <- list(c("control cortex","control WM"),c("control cortex","myelinated cortex"),
                       c("control cortex","demyelinated cortex"),c("control WM","myelinated cortex"),
                       c("control WM","demyelinated cortex"),c("myelinated cortex","demyelinated cortex"))
colors <- brewer.pal(n=4,name = "Dark2")
# colors <- colors[c(1,3,2)]
gene <- "SDC2"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "pathology", cols = colors,
              # idents = c("Microglia")
              idents = c("16","18","34","38")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.00001, max(p1[[1]][["data"]][[gene]])+0.8*max(p1[[1]][["data"]][[gene]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  xlab("")
p1$layers[[2]]$aes_params$alpha <- 0.3
p1



eicosanoids <- c("PTGS2","FADS1","FADS2","FADS3","EGR2","EGR3","PLA2G10","PLA2G4A","DAGLB","PLA2G4E",
                 "PLA2G3", "PLBD1","ACSL1","ACSL2","ACSL3","ACSL4",
                 "PLA2G4F","PLA2G6","PLA2G7","PLA2G5","PLA2G4C","ABCC4","PTGES2","SLCO2A1")

plcs <- c("ABCC4","SLCO2A1",
          "PTGS1","PTGS2","PTGES2",
          "PTGES","PTGES3",
          "ALOX5","ALOX5AP","ALOX12","ALOX15",
          "ELOVL2","FADS1","ELOVL5","FADS2","FADS3",
          "DAGLA","DAGLB",
          "PLA2G5","PLA2G2F","PLA2G12A","PLA2G10","PLA2G12B","PLA2G4A","PLA2G4E","PLA2G7","PLA2G4C","PLA2G6",
          "PLA2G2D","PLA2G2E","PLA2G2C","PLA2G3","PLA2G1B","PLA2G2A","PLA2G4F","PLBD1",
          "ACSL3","ACSL4","ACSL1", "ACSL2"
)

plcs2 <- c("ACSL3","ACSL4","ACSL1", # "ACAT1","ACAT2",
           # "PLA2G7",
           "FADS3","PTGS2", # "FADS1","FADS2","ALOX5AP","LTA4H",
           "PTGES3","TBXAS1"
)

ptgs <- c("ABCC4","SLCO2A1","PTGES2","PTGS2")

eicosanoids2 <- c("ABCC4","PTGES2","PTGS2","FADS1","ACSL1","PLA2G10","PLA2G4A","PLBD1")

genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")
clus4.markers <- c("SDC2","PTGS2","HBEGF","IL1B","FOSL2","EGR2","EGR3","CXCL8", # "PPIF","NLRP3","NFKB1","GPR183","EGR1","RGCC","TRIB1",
                   "G0S2")


DotPlot(prot.combined, features = plcs2,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        idents = c("16","18","34","38"),
        cluster.idents = F, 
        group.by = "pathology",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + xlab("") + ylab("")




############### DEG analysis

cell.types <- c(
  "Microglia"
)
# cell.types <- as.character(1:7)

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
prot.combined$celltypes.disease <- paste(prot.combined$type_broad, prot.combined$disease, sep = "_")
prot.combined$celltypes.lesion <- paste(prot.combined$type_broad, prot.combined$lesion_type, sep = "_")
prot.combined$cluster.lesion <- paste(prot.combined$seurat_clusters, prot.combined$lesion_type, sep = "_")
# Idents(prot.combined) <- "cluster.disease"
Idents(prot.combined) <- "celltypes.lesion"


# ## wilcox
lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             # ,"CAL" #,"CIL","RL",
)
up.list <- c()
down.list <- c()

for (i in 1:length(cell.types)) {
  for ( e in 1:length(lesions)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_CIL"),
                                ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                slot = "data",
                                assay = "RNA",
                                features = NULL,
                                logfc.threshold = 0,
                                test.use = "wilcox",
                                min.pct = 0.0,
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
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_CIL_x_",lesions[e],"_", cell.types[i], "_by.age.sex_DEGs.xlsx"))
    zk.response1 <- zk.response0[zk.response0$avg_log2FC > 0.322,]
    zk.response2 <- zk.response0[zk.response0$avg_log2FC < -0.322,]
    gc()
    zk.response0b <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_CAL"),
                                 ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                 slot = "data",
                                 assay = "RNA",
                                 features = NULL,
                                 logfc.threshold = 0,
                                 test.use = "wilcox",
                                 min.pct = 0.0,
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
    zk.response0b <- zk.response0b[zk.response0b$p_val_adj < 0.1,]
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_CAL_x_",lesions[e],"_", cell.types[i], "_DEGs.xlsx"))
    zk.response1b <- zk.response0b[zk.response0b$avg_log2FC > 0.322,]
    zk.response2b <- zk.response0b[zk.response0b$avg_log2FC < -0.322,]
    gc()
    zk.response0c <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_RL"),
                                 ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                 slot = "data",
                                 assay = "RNA",
                                 features = NULL,
                                 logfc.threshold = 0,
                                 test.use = "wilcox",
                                 min.pct = 0.0,
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
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_RL_x_",lesions[e],"_", cell.types[i], "_DEGs.xlsx"))
    zk.response1c <- zk.response0c[zk.response0c$avg_log2FC > 0.322,]
    zk.response2c <- zk.response0c[zk.response0c$avg_log2FC < -0.322,]
    gc()
    zk.response1e <- zk.response1b[row.names(zk.response1b) %in% row.names(zk.response1),]
    zk.response1e <- zk.response1e[row.names(zk.response1e) %in% row.names(zk.response1c),]
    zk.response2e <- zk.response2b[row.names(zk.response2b) %in% row.names(zk.response2),]
    zk.response2e <- zk.response2e[row.names(zk.response2e) %in% row.names(zk.response2c),]
    up.list[[e]] <- row.names(zk.response1e)
    down.list[[e]] <- row.names(zk.response2e)
    rm(zk.response1)
    rm(zk.response2)
    rm(zk.response0)
    rm(zk.response1b)
    rm(zk.response2b)
    rm(zk.response0b)
    rm(zk.response1c)
    rm(zk.response2c)
    rm(zk.response0c)
    rm(zk.response1e)
    rm(zk.response2e)
    gc()
  }
}

names(up.list) <- lesions
up.list2 <- t(plyr::ldply(up.list, rbind))
colnames(up.list2) <- up.list2[1,]
up.list2 <- up.list2[-c(1), ] 
write.xlsx(as.data.frame(up.list2), rowNames = F,file="wilcox_microglia_chronic_lesions_converged_up_per.non.chronic_DEGs.xlsx")
names(down.list) <- lesions
down.list2 <- t(plyr::ldply(down.list, rbind))
colnames(down.list2) <- down.list2[1,]
down.list2 <- down.list2[-c(1), ] 
write.xlsx(as.data.frame(down.list2), rowNames = F,file="wilcox_microglia_chronic_lesions_converged_down_per.non.chronic_DEGs.xlsx")

lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             # ,"CAL" #,"CIL","RL",
)
common_up_genes <- Reduce(intersect, list(up.list2[,1],up.list2[,3],up.list2[,4],up.list2[,5],up.list2[,6]))



for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = c(paste0(cell.types[i], "_CAL"),
                                                    paste0(cell.types[i], "_CIL"),
                                                    paste0(cell.types[i], "_RL")),
                              ident.2 = c(paste0(cell.types[i], "_AL"),
                                          paste0(cell.types[i], "_GM"),
                                          paste0(cell.types[i], "_GML"),
                                          paste0(cell.types[i], "_NAGM"),
                                          paste0(cell.types[i], "_NAWM"),
                                          paste0(cell.types[i], "_WM")),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "wilcox",
                              min.pct = 0.0,
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_chronic_lesions_x_all_others_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}


## MAST

lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             #,"CAL" #,"CIL","RL",
)
for (i in 1:length(cell.types)) {
  for ( e in 1:length(lesions)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_RL"),
                                ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                slot = "data",
                                assay = "RNA",
                                features = NULL,
                                logfc.threshold = 0,
                                test.use = "MAST",
                                min.pct = 0.0,
                                min.diff.pct = -Inf,
                                verbose = TRUE,
                                only.pos = FALSE,
                                max.cells.per.ident = Inf,
                                random.seed = 1,
                                latent.vars = c("sex","age_cat"),
                                min.cells.feature = 3,
                                min.cells.group = 3,
                                pseudocount.use = 0.001,
                                mean.fxn = NULL,
                                fc.name = NULL,
                                base = 2,
                                densify = FALSE,
                                #recorrect_umi = TRUE
    )
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_RL_x_",lesions[e],"_", cell.types[i], "_by.age.sex_DEGs.xlsx"))
    rm(zk.response0)
    gc()
  }
}


for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_CAL"),
                              ident.2 = paste0(cell.types[i], "_AL"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("sex","age_cat"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_CAL_x_AL_", cell.types[i], "_minpct0_DEGs.xlsx"))
  rm(zk.response0)
  gc()
  # zk.response0 <- FindMarkers(myeloids, ident.1 = c(paste0(cell.types[i], "_CAL"),paste0(cell.types[i], "_CIL"),paste0(cell.types[i], "_RL")),
  #                             ident.2 = paste0(cell.types[i], "_HC"),
  #                             slot = "data",
  #                             assay = "RNA",
  #                             features = NULL,
  #                             logfc.threshold = 0,
  #                             test.use = "MAST",
  #                             min.pct = 0.0,
  #                             min.diff.pct = -Inf,
  #                             verbose = TRUE,
  #                             only.pos = FALSE,
  #                             max.cells.per.ident = Inf,
  #                             random.seed = 1,
  #                             latent.vars = c("sex","age_cat"),
  #                             min.cells.feature = 3,
  #                             min.cells.group = 3,
  #                             pseudocount.use = 0.001,
  #                             mean.fxn = NULL,
  #                             fc.name = NULL,
  #                             base = 2,
  #                             densify = FALSE,
  #                             #recorrect_umi = TRUE
  # )
  # zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_PMS_x_HC_clus", cell.types[i], "_minpct0.1_DEGs.xlsx"))
  # rm(zk.response0)
  # gc()
  # zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_RRMS"),
  #                             ident.2 = paste0(cell.types[i], "_HC"),
  #                             slot = "data",
  #                             assay = "RNA",
  #                             features = NULL,
  #                             logfc.threshold = 0,
  #                             test.use = "MAST",
  #                             min.pct = 0.1,
  #                             min.diff.pct = -Inf,
  #                             verbose = TRUE,
  #                             only.pos = FALSE,
  #                             max.cells.per.ident = Inf,
  #                             random.seed = 1,
  #                             latent.vars = c("sex"),
  #                             min.cells.feature = 3,
  #                             min.cells.group = 3,
  #                             pseudocount.use = 0.001,
  #                             mean.fxn = NULL,
  #                             fc.name = NULL,
  #                             base = 2,
  #                             densify = FALSE,
  #                             #recorrect_umi = TRUE
  # )
  # zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_RRMS_x_HC_clus", cell.types[i], "_minpct0.1_DEGs.xlsx"))
  # rm(zk.response0)
  # gc()
}


########## subset microglia

myeloids <- subset(prot.combined, subset = celltypes %in% c("Microglia"))

# myeloids <- subset(prot.combined, subset = unintegrated_clusters %in% c("24","29","26","31","40","46","59","71"))

### Normalize and scale merged obj
myeloids <- NormalizeData(myeloids, normalization.method = "LogNormalize")
myeloids <- FindVariableFeatures(myeloids, selection.method = "vst", nfeatures = 1000)
myeloids <- ScaleData(myeloids, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                      model.use = "linear")#, features = all.genes)
gc()

### Dimensionality reduction and integration
myeloids <- RunPCA(myeloids, npcs = 50)
gc()
ElbowPlot(myeloids, ndims = 50)
myeloids <- FindNeighbors(myeloids, dims = 1:30, reduction = "pca")
myeloids <- FindClusters(myeloids, resolution = 0.1, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
myeloids <- RunUMAP(myeloids, #umap.method = "umap-learn", 
                    dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(myeloids, reduction = "umap.unintegrated", raster = F,
        ncol = 2,
        label = F,
        repel = T,
        #group.by = "patient"
        split.by = "pathology"
)#, combine = F)

myeloids <- BuildClusterTree(myeloids, assay = "RNA", reduction = "umap.unintegrated", reorder = T)
prot.markers <- FindAllMarkers(myeloids, assay = "RNA", only.pos = T, min.pct = 0.2, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.clusters.myeloids.markers.xlsx")

DotPlot(myeloids, features = c("P2RY12","CSF1R","CD74","C3","CST3","C1QA", "CX3CR1","TMEM119","SLC2A5",
                               "IL1B","NFKB1","APOE","IRF8", "CST7","CLEC7A","AIF1",# "HEXB", 
                               "MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15",
                               "EMILIN2", "MARCO", #"CD4","GDA", "PLA2G7","PLXND1",#"ASAH1","GAS7",
                               "PPIF","PTGS2","PTGES3","NLRP3","HBEGF","FOSL2","GPR183",#"TBXAS1",
                               "CXCL8","G0S2","RGCC","SDC2" #,"TRIB1" #"EGR1","EGR2","EGR3",
                               # "CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7"#,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
),
col.max = 20,
dot.scale = 10, 
cluster.idents = T, #group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

myeloids <- RenameIdents(myeloids, 
                         `0`= "Homeostatic microglia",
                         `1`= "Foam cells",
                         `2`= "Monocytes",
                         `3`= "Inflammatory microglia",
                         `4`= "Homeostatic microglia",
                         `5`= "Homeostatic microglia",
                         `6`= "Homeostatic microglia",
                         `7`= "Homeostatic microglia"
)
myeloids$celltypes <- sub("(.*)_.*", "\\1", Idents(myeloids))
Idents(myeloids) <- "celltypes"


## Save and Load data
saveRDS(object = myeloids, file = "myeloids_reg.out.mito.ncounts.patients.annotated.v3.Rds")
# #rm(myeloids)
myeloids <- readRDS("./myeloids_reg.out.mito.ncounts.patients.annotated.v3.Rds")

##### freq
num.cells <- as.data.frame(table(myeloids$disease_lesion, myeloids$disease_lesion))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$disease_lesion, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.myeloids_disease_lesion.csv")

paletteLength <- 40
colors <- colorRampPalette( rev(brewer.pal(11, "Set3")))(paletteLength)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
barplot(tfreq.num.cells.celltype, col = sample(col_vector, n), legend.text = rownames(tfreq.num.cells.celltype),
        xlim = c(0,15.5), main = "Myeloid cell types distribution", xlab = "", ylab = "Cells percentage (%)", cex.names = 0.9)

myeloids$lesion.patient <- paste(myeloids$lesion_type, myeloids$GEM, sep = "_")
num.cells <- as.data.frame(table(myeloids$lesion.patient, myeloids$lesion.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$lesion.patient, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.myeloids_lesion_patient.csv")

mydata <- read.csv("freq.num.cells.myeloids_lesion_patient.csv", header = T)

coluna1 <- c("GM","GML","NAGM",
             "WM","AL","NAWM",
             "CAL","CIL","RL")

my_comparisons <- list(#c("RL","GM"),c("RL","GML"),c("RL","NAGM"),c("RL","WM"),c("RL","AL")
  c("CAL","GM"),c("CAL","GML"),c("CAL","NAGM"),c("CAL","WM"),c("CAL","AL")
  # c("CIL","GM"),c("CIL","GML"),c("CIL","NAGM"),c("CIL","WM"),c("CIL","AL")
  # c("RL","GM"),c("RL","GML"),c("RL","NAGM"),c("RL","WM"),c("RL","AL")
)

p1 <- ggplot(mydata, aes_string(x="factor(lesion,levels=coluna1)", y="Foam.cells", color="lesion")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Set1") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("cells percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  # ggtitle("COX2+ Classical monocytes") +
  # ggtitle("Nonclassical monocytes") +
  ggtitle("Foam cells") +
  # ggtitle("Intermediate monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1


clus4.markers <- c("SDC2","PPIF","PTGS2","NLRP3","HBEGF","IL1B","NFKB1","FOSL2","GPR183","EGR1","EGR2","EGR3","CXCL8",
                   "G0S2","RGCC","TRIB1")

clus4.markers <- c("PPIF","PTGS2","PTGES3","NLRP3","HBEGF","FOSL2","GPR183",
                   "CXCL8","G0S2","RGCC","SDC2")

myeloids <- AddModuleScore(
  myeloids,
  features = clus4.markers,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "clus4.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)


gene <- "SDC2"
FeaturePlot(myeloids, features = gene, 
            raster = T,
            pt.size = 2.,
            # ncol = 3,
            reduction = "umap.unintegrated", 
            # min.cutoff = 0.7,
            max.cutoff = 2,
            cols = c("gray90","#9E1021"),
            split.by = "pathology"
) + theme(legend.position = "right")

# Idents(myeloids) <- "celltypes"
Idents(myeloids) <- "type_broad"

# my_comparisons <- list(c("PMS","RRMS"),c("PMS","HC"),c("RRMS","HC"))
# colors <- brewer.pal(n=9,name = "Dark2")
# colors <- colors[c(1,3,2)]
gene <- "PTGES3"
p1 <- VlnPlot(myeloids, features = gene,
              pt.size = 0.05, raster = F, 
              group.by = "lesion_type", #cols = colors,
              # group.by = "disease_lesion", 
              # idents = c("Microglia")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 20, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000000, max(p1[[1]][["data"]][[gene]])+0.*max(p1[[1]][["data"]][[gene]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1



eicosanoids <- c("PTGS2","FADS1","FADS2","FADS3","EGR2","EGR3","PLA2G10","PLA2G4A","DAGLB","PLA2G4E",
                 "PLA2G3", "PLBD1","ACSL1","ACSL2","ACSL3","ACSL4",
                 "PLA2G4F","PLA2G6","PLA2G7","PLA2G5","PLA2G4C","ABCC4","PTGES2","SLCO2A1")

plcs <- c("ABCC4","SLCO2A1",
          "PTGS1","PTGS2","PTGES2",
          "PTGES","PTGES3",
          "ALOX5","ALOX5AP","ALOX12","ALOX15",
          "ELOVL2","FADS1","ELOVL5","FADS2","FADS3",
          "DAGLA","DAGLB",
          "PLA2G5","PLA2G2F","PLA2G12A","PLA2G10","PLA2G12B","PLA2G4A","PLA2G4E","PLA2G7","PLA2G4C","PLA2G6",
          "PLA2G2D","PLA2G2E","PLA2G2C","PLA2G3","PLA2G1B","PLA2G2A","PLA2G4F","PLBD1",
          "ACSL3","ACSL4","ACSL1", "ACSL2"
)

plcs2 <- c("ACSL3","ACSL4","ACSL1",
           "PLA2G7",
           "FADS3", #"FADS2","FADS1",
           "ALOX5AP","LTA4H","PTGES3","TBXAS1"
)

ptgs <- c("ABCC4","SLCO2A1","PTGES2","PTGS2")

eicosanoids2 <- c("ABCC4","PTGES2","PTGS2","FADS1","ACSL1","PLA2G10","PLA2G4A","PLBD1")

genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")


DotPlot(myeloids, features = plcs2,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        # idents = "Microglia",
        cluster.idents = F, 
        group.by = "lesion_type",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


















