library(Seurat)
library(scCustomize)
library(ggplot2)
library(SCP)
# BiocManager::install("dittoSeq")
library(dittoSeq)
library(ggpubr)


setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
# setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
gc()
# merged_females <- readRDS("data/merged_seurat/females/females_merged_object.rds")

merged <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")

# genes for violin plots
# FOS, FOSL2, PTGS2, PTGES3, TBXAS1

# for mo-DCs
# CD83, PTGER4, IL1b, RIPK2, NAMPT, BTG2, ATF, NFKBIZ, PLAUR


DS_cols <- c("#7ba39d", "#cab2d6","#fdbf6f", "#a6cee3")

vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(merged, features = signature,
            pt.size = 0.1, 
            group.by = "prog", 
            y.max = y_max,
            cols = DS_cols# add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(method = "anova", label.y = 40)+ # Add global p-value
      stat_compare_means(aes(label = after_stat(p.signif)),
                         method = "t.test", ref.group = "0.5")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 12, height = 18)
}

gene_sig <- c("FOS", "FOSL2", "PTGS2", "PTGES3", "TBXAS1")
comparisons <- list(c("HC", "non-progressor",
                      "progressor", "RR"))
vp_case1(gene_signature = gene_sig, file_name = "figures_paper2/fig3/", test_sign = comparisons)

table(merged$prog)

?stat_compare_means


markers_test <- FindMarkers(merged, features = gene_sig, group.by = "prog", ident.1 = "HC")

?FindMarkers
