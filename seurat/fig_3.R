library(Seurat)
library(scCustomize)
library(ggplot2)
library(SCP)
# BiocManager::install("dittoSeq")
library(dittoSeq)


setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
# setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
gc()
# merged_females <- readRDS("data/merged_seurat/females/females_merged_object.rds")

merged <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")

clust_colors <- c(
  "#E6B8B8",  # Darkened #F4CCCC
  "#CFCFEF",  # Darkened #E6E6FA
  "#EFE49A",  # Darkened #FFFACD
  "#C185C1",  # Darkened #DDA0DD
  "#8EB5D6",  # Darkened #ADD8E6
  "#9ACB6F",  # Darkened #b2df8a
  "#E6C4C9",  # Darkened #FFF0F5
  "#B0B0B0",  # Darkened #D3D3D3
  "#D4E1F1",  # Darkened #F0F8FF
  "#E6C099",  # Darkened #FFDAB9
  "#6EA8CC",  # Darkened #87CEEB
  "#E6C092",  # Darkened #FFE5B4
  "#A0CFA1"   # Darkened #C7E5C9
)


table(merged$prog)

Idents(merged) <- "sex"
males <- subset(merged, idents = "M", invert = T)
clust_split_males <- DimPlot(males, group.by = "seurat_clusters", split.by = "prog", raster = F,
                       cols = clust_colors, label = T)
clust_split_males

ggsave("figures_paper2/fig3/females_clust_UMAP_split.pdf", device = "pdf", 
       width = 18, height = 8, units = "in", dpi = 100, plot = clust_split_males)
