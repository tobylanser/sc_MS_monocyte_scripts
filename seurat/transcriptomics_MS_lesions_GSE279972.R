#############################################
####################
########## starting from raw abundance

############## transcriptomics


# libraries
library(grid)
library(ggstatsplot)
library(data.table)
library(ggpubr)
library(car)
library(pheatmap)
library(openxlsx)
library(RColorBrewer)
library(ggplot2)
library(matrixStats)
library(tidyverse)
library(zoo)
library(EnhancedVolcano)
library(missForest)     # RF imputation 
library(imputeLCMD)     # QRILC, KNN, MAR/MNAR helpers
library(BiocGenerics)
library(tximport)
library(S4Vectors)
library(DESeq2)
library(biomaRt)
library(apeglm)
library(readr)
library(robustbase)
library(corrplot)
library(PerformanceAnalytics)
library(genefilter)
library(rnaseqGene)
# library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(BiocParallel)
register(MulticoreParam(14))



setwd("/run/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/bulkRNAseq_MS_lesions_GSE279972/transcriptomics/")

data1 <- read.xlsx("Processed_data_all_omics.xlsx", sheet = 4)


### make transcript to gene file

# mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")#, host="https://asia.ensembl.org")
# mart <- useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl", #mirror = "asia")
#                    host="https://www.ensembl.org")
# 
# 
# txinfo <- getBM(attributes = c('ensembl_transcript_id_version',
#                                'ensembl_transcript_id', 
#                                'ensembl_gene_id', 
#                                'external_transcript_name',
#                                'external_gene_name'),
#                 #filters = 'ensembl_transcript_id', 
#                 #values = transcript_ids,
#                 mart = mart)%>%
#   as_tibble() %>%
#   mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))
# 
# tx2gene0 <- txinfo %>% dplyr::select(external_gene_name, ensembl_gene_id) %>%
#   dplyr::rename(c("Gene0" = "external_gene_name", "Gene_ID" = "ensembl_gene_id" )) %>%
#   data.frame()

txinfo <- transcripts(EnsDb.Hsapiens.v86, 
                      columns = c("tx_id", "gene_id", "tx_name", "gene_name")) %>%
  as_tibble()

tx2gene0 <- txinfo %>% dplyr::select(gene_name, gene_id) %>%
    dplyr::rename(c("Gene0" = "gene_name", "Gene_ID" = "gene_id" )) %>%
    data.frame()
tx2gene0$Gene <- with(tx2gene0, ifelse(Gene0 == "", Gene_ID, Gene0))
tx2gene <- tx2gene0[!duplicated(tx2gene0$Gene_ID), ]
rm(tx2gene0)
tx2geneB <- tx2gene[tx2gene$Gene_ID %in% data1$X1,]

data1$Gene <- data1$X1

for (i in 1:length(tx2geneB$Gene_ID)) {
  data1$Gene <- car::recode(data1$Gene, "tx2geneB$Gene_ID[i] = tx2geneB$Gene[i]")
}

data1_unique <- data1[!duplicated(data1$Gene), ]
row.names(data1_unique) <- data1_unique$Gene
mydf_Mprime <- data1_unique[,2:106]

meta <- read.xlsx("Processed_data_all_omics.xlsx", sheet = 2)
meta$ID <- gsub(" ",".", meta$Unifying_code)
# meta$ID <- gsub("-",".", meta$ID0)
sampleTable <- data.frame(sampleName = meta$ID,
                          sex = meta$sex,
                          age = meta$age,
                          morphology_mic = meta$Morphology.microglia,
                          Lesion_type_8 = meta$Lesion_type_8,
                          Lesion_type_9 = meta$Lesion_type_9,
                          Lesion_type_6 = meta$Lesion_type_6,
                          disease = meta$diagnosis,
                          donor = meta$NBB.donor.ID
)

## Ensure basic structure
stopifnot(all(colnames(mydf_Mprime) %in% sampleTable$sampleName))

## Align sampleTable to expression columns
sampleTable_sub <- sampleTable[match(colnames(mydf_Mprime), sampleTable$sampleName), ]
stopifnot(all(sampleTable_sub$sampleName == colnames(mydf_Mprime)))

## Factors: choose reference levels so that the ProgMS–StRRMS contrast is explicit
sampleTable_sub$disease  <- factor(sampleTable_sub$disease, levels = c("Non-demented control","Multiple sclerosis"))
sampleTable_sub$sex <- factor(sampleTable_sub$sex, levels = c("f", "m"))
sampleTable_sub$age <- as.numeric(sampleTable_sub$age)  # continuous covariate
sampleTable_sub$Lesion_type_6 <- factor(sampleTable_sub$Lesion_type_6, levels = c("CWM","NAWM","2","3","4","6"))
sampleTable_sub$Lesion_type_9 <- factor(sampleTable_sub$Lesion_type_9, levels = c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2",
                                                                                  "3.3","4","6"))

## Container for results
res <- data.frame(
  accession = rownames(mydf_Mprime),
  diff_exp = NA_real_,
  p_val = NA_real_,
  diffexp_sexm = NA_real_,
  p_val_sexm = NA_real_,
  stringsAsFactors = FALSE
)

group_compare <- "4"

### run ANCOVA
## Loop over metabolites, fit lm: expression ~ group + gender
for (i in seq_len(nrow(mydf_Mprime))) {
  y <- as.numeric(mydf_Mprime[i, ])
  ## Skip metabolites with all-NA or too many NAs
  if (all(is.na(y))) next
  df <- data.frame(
    expr   = y,
    lesion  = sampleTable_sub$Lesion_type_6,
    sex = sampleTable_sub$sex,
    age = sampleTable_sub$age
  )
  ## Remove samples with missing expression
  # df <- df[!is.na(df$expr), ]
  df <- df[complete.cases(df), ]
  if (nrow(df) < 5) next  # too few samples to fit
  fit <- try(lm(expr ~ lesion + sex + age, data = df), silent = TRUE)
  # fit <- try(lm(expr ~ group + gender, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) next
  coefs <- summary(fit)$coefficients
  ## With StRRMS as reference, the ProgMS–StRRMS difference is the coefficient of groupProgMS
  if (paste0("lesion",group_compare) %in% rownames(coefs)) {
    res$diff_exp[i]  <- coefs[paste0("lesion",group_compare), "Estimate"]
    res$p_val[i] <- coefs[paste0("lesion",group_compare), "Pr(>|t|)"]
  }
  res$diffexp_sexm[i]  <- coefs["sexm", "Estimate"]
  res$p_val_sexm[i] <- coefs["sexm", "Pr(>|t|)"]
}


########

my_res <- res

# Export results
write.xlsx(as.data.frame(my_res), rowNames = T, file="RNAseq_Lesion6_IL_vs_CWM_stats_lowess_sex_age_ANCOVA_imputation_stats.xlsx")


## PCA for QC
# rld <- lg2mydfmean
rld <- mydf_Mprime
# rld <- mydf_Mprime[,c(1,2,5,6,11,12,15,16,23,26,3,4,9,10,17,18,21,22,24,25)]
ntop <- 200
rv <- rowVars(as.matrix(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(rld[select, ]))
sample_names <- colnames(rld)
# sample_names <- my_samp
# sample_names <- my_samp[c(1,2,5,6,11,12,15,16,23,26,3,4,9,10,17,18,21,22,24,25)]
group <- sampleTable_sub$Lesion_type_6
# group <- factor(group)
group <- factor(group, levels = c("CWM","NAWM","2","3","4","6"))#,
# group <- factor(sampleTable$Lesion_type_9, levels = c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2",
#                                                       "3.3","4","6"))#,
# "ProgMS_G2","ProgMS_GL2","ProgMS_PM"))
# group <- factor(sampleTable$disease, levels = c("CTRL","ProgMS"))
# group <- factor(ifelse(my_condition == "RRMS", "RRMS", "ProgMS"))
# sample_names <- colnames(mydf_Mprime)
# group <- c("RRMS","RRMS","RRMS","ProgMS","ProgMS","ProgMS","ProgMS","ProgMS")
# group <- factor(group, levels = c("RRMS","ProgMS"))
plot(pca$x[,1], pca$x[,2],
     col=group, pch=19, xlab="PC1", ylab="PC2",
     main="PCA of Samples")
legend("bottomright", legend=levels(group), col=1:length(levels(group)), pch=19)
text(pca$x[,1], pca$x[,2], labels=sample_names, pos=2, cex=0.7)

### heatmap DEP

# Identify DEPs
dep_1 <- my_res[!is.na(my_res$p_val),]
dep_2 <- dep_1[dep_1$p_val < 0.05 , ]
# dep_2 <- dep_2[dep_2$diff_exp > 0 , ]
# dep_2 <- dep_2[abs(dep_2$diff_exp) < 0.3 , ]
# dep_2 <- dep_1[dep_1$stat == TRUE & abs(dep_1$diffexp) > 0.5, ]
# dep_2 <- dep_1[dep_1$stat == TRUE , ]
# dep_2 <- dep_2[!duplicated(dep_2$Gene), ]
# heatmap_data <- lg2mydfmean[row.names(lg2mydfmean) %in% row.names(dep_2), ]
plcs <- c(
  #"SLCO2A1","PTGS1","PTGES",
  "ABCC4","PTGES3","PTGES2","PTGS2",
  #"ALOX5","ALOX5AP","ALOX12","ALOX15",
  "ELOVL5","FADS1","FADS2","FADS3","ELOVL2",
  #"DAGLA","DAGLB","PLA2G5","PLA2G2F","PLA2G12A","PLA2G10","PLA2G12B","PLA2G4A","PLA2G4E", 
  #"PLA2G4C","PLA2G6","PLA2G2D","PLA2G2E","PLA2G2C","PLA2G3","PLA2G1B","PLA2G2A","PLA2G4F",
  "PLBD1","PLA2G7",
  "ACSL3","ACSL1" #, "ACSL4","ACSL2","MBOAT7" 
)
genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")

carnitine <- c("SLC22A5","CPT2","CPT1A","SLC25A20", # carnitine
               "ADIPOR1","PLIN2","PPT1","ASAH1", # lipid metabolism and storage
               "PLSCR1", # phosphatidyl serine
               "CBX3","GPX4","TFRC","SLC25A37","STEAP4","SLC11A1",  # ferroptosis
               "UGCG",  # glycosphingolipids
               "PGS1", # phosphatidyl glycerol
               "TSPO","TSPOAP1","STARD4","SREBF2","LDLR","HMGCR","APOM","SREBF1","ABCA1","ABCG1","ABCG5","ABCG8","NPC1", # cholesterol
               "LIPA","CES1",  # cholestery ester
               "FLVCR2","LPCAT1","LPCAT2","LPCAT3","LPCAT4","GPCPD1","CHKA","GPAT4"  # choline
)  
# TG <- c("","","","","")
choline  <- c("SCARB1","MFSD2A","FABP3","LPCAT1","LDLR","APOC1","RAB38","ACSL3","ABCA3","CAPN2")
lipids <- c("ACSL1","ACSL3","ACSL5","PLIN2")

genes1 <- c("CPT2","CPT1A","SLC25A20","LCLAT1","HADHA","ALDH9A1", # carnitine
            "ADIPOR1", # lipid metabolism and storage
            "PLSCR1", # phosphatidyl serine
            "GPX4","TFRC","SLC25A37","STEAP4","SLC11A1","SLC7A11",  # ferroptosis
            # "UGCG",  # glycosphingolipids
            "PGS1", # phosphatidyl glycerol
            # "TSPO","TSPOAP1",#"STARD4",#"SREBF2", # cholesterol
            "FLVCR2","LPCAT1","CHKA","CHKB","ETNK1","ETNK2","PCYT1A","PCYT1B","PCYT2","CEPT1","SELENOI","GPAT4","PLA2G4A","PTDSS1"  # choline "GPCPD1",,"LCAT","PLB1","PLA2G7","PLCE1","PLCH1","ATP8A2"
)  
lands <- c("MBOAT1","MBOAT2","HHATL","MBOAT4","MBOAT7","LPCAT4","LPCAT3","LPCAT2","LPCAT1","PLA2G4A") ## lands cycle

phospholipases <- c("DAGLA","DAGLB","PLA2G5","PLA2G2F","PLA2G12A","PLA2G10","PLA2G12B","PLA2G4A","PLA2G4E", 
            "PLA2G4C","PLA2G6","PLA2G2D","PLA2G2E","PLA2G2C","PLA2G3","PLA2G1B","PLA2G2A","PLA2G4F",
            "PLBD1","PLA2G7") # phospholipases

kennedy <- c("CHPT1","CHKB","ETNK1","ETNK2","PCYT1A","PCYT1B","PCYT2","CEPT1","SELENOI") ## kennedy pathway

genes1 <- c("SLC25A48","CHDH","CHAT","SLC18A3") # "FLVCR2","ALDH9A1" 

genes1 <- c("MFSD2A","MFSD2B","SPNS1","SPNS2","ATP8B1","SGMS1","SGMS2","SMPD1")  ## sphingomyelin metabolism

peroxisome <- c("CYBA","PARK7","RAC2","ECH1","GSTK1","DHRS4","RPS27A", ## peroxisome related associated with PIRA2  "UBE2D1","RAC1","ACSL1","ACAA1","RHOB","PEX16","ACSL3",
            "UBA52","UBB","FPR2","GNAI2","PMVK", 
            "GSTP1","ITGAM","ITGB2","TYROBP","CLEC7A","CFL1","CYP1B1","DUSP1",
            "ERN1","RPS3","TXN","KLF4","PPIF","KDM6B","SIRPA","PRDX5","PRDX4","PRDX3","PRDX2","PRDX1","SOD1"
) # "SOD2","SCP2","HSD17B4","CYBB","IDH1","UBE2D3","FYN","TRPM2","OSER1","EEF2","KLF2",

peroxisome2 <- c("PARK7","RAC2","RPS27A","UBA52","FPR2", ## peroxisome related associated with foam cells cluster "ACSL1","ACAA1","RHOB","ACSL3",
            ## "SOD1","PRDX5","PRDX4","PRDX3","RPS3","ITGAM","ITGB2","TYROBP","PEX16","PMVK","GNAI2","UBB","UBE2D1","CYBA","ECH1","GSTK1","DHRS4","RAC1",
            "GSTP1","CLEC7A","CFL1","CYP1B1","DUSP1","ERN1","TXN","KLF4","PPIF","KDM6B","SIRPA","PRDX2","PRDX1"
) # "SOD2","SCP2","HSD17B4","CYBB","IDH1","UBE2D3","FYN","TRPM2","OSER1","EEF2","KLF2",

polyamines <- c("AMD1","OAZ1","ODC1","SAT1","SRM","SMS","SMOX","DHPS","DOHH","EIF5A")  ## polyamines

phagocytosis <- c("LYST","NCF2","BIN2","RAB14","CLCN3","PRKCD","SPG11","ITGAL","P2RX7","PECAM1", ## phagocytosis "CD93","CD36","IRF8",
            "TICAM2","VAV1","CCR2","ICAM3","PAK1","ARHGAP25","ITGB2","ABL1", ### phagocytosis  "FCN1","ANXA1","ITGAM","CORO1A","RAC1","FCGR1A","NCF4","CD14","VAV2","HCK",
            "ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0A2",  ## ROS "ATP6V0A1",
            "LIPA","DNAJC13","HEATR5A","DENND1A","FKBP15","LRRK2","RIN2","SNX10"  ## endocytosis  "ASGR1","PYCARD","SNX17","CAP1",
)

endocytosis <- c("CD93","LYST","SNX17","CLCN3","CORO1A","ITGAL","P2RX7","PECAM1",
                 "LIPA","DNAJC13","FCGR1A","CAP1","CCR2","HEATR5A","ASGR1","PYCARD",
                 "DENND1A","FKBP15","ITGB2","HCK","LRRK2","CD36","RIN2","SNX10")
ROS <- c(
  # "CYBB","RAC2","HVCN1","NCF1",
  "ATP6V1D","ATP6V1F","ATP6V1G1","ATP6V1G3","ATP6V1C2","ATP6V1B1","ATP6V1E2","ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1H",
  "ATP6V0A3","ATP6V0C","ATP6V0E1","ATP6V0E2","ATP6V0D2","ATP6V0A4" #,"ATP6V0A1","ATP6V0A2","ATP6V0D1","ATP6V1G2","ATP6V1E1",
)

ferroptosis <- c("FLVCR2","TFRC","FTL","SLC39A8","PCBP2", "NCOA4","ATG7","ATG5","SLC40A1",#"FTH1",
            "SLC25A37","STEAP4","STEAP3","SLC11A1","SLC7A11","TMEM164","GSS","GCLM","GPX4"  # ferroptosis # "SLC33A1",
)
clus4.markers <- c("SDC2","SDC4","G0S2","EGR1","HBEGF","FOSL2","IL1B","PTGS2","EGR2","EGR3","CCL3L3","IRAK2","DUSP2","B3GNT5","TSPOAP1",
                   "CXCL8","TRIB1","PPIF","NLRP3","GPR183" #, "NFKB1","RGCC",
)
clus4.markers <- unique(c("SDC2","SDC4","G0S2","EGR1","FOSL2","IL1B","PTGS2","EGR2","EGR3","CCL3L3","DUSP2","B3GNT5","TSPOAP1","CXCL8",
                   "TRIB1","NLRP3","GPR183","LPCAT1",ROS,peroxisome,peroxisome2,ferroptosis # "PPIF","HBEGF","IRAK2","NFKB1","RGCC","KLF15",
))

# targets0 <- unique(c(ferroptosis,lands,ROS,endocytosis,phagocytosis,phospholipases))

targets <- clus4.markers[clus4.markers %in% dep_2$accession]
# targets <- grep("^PC\\(.*", dep_2$accession, value = T)
# targets <- grep("^PC\\(O-.*", dep_2$accession, value = T)
heatmap_data <- mydf_Mprime[row.names(mydf_Mprime) %in% targets, ]
# heatmap_data <- mydf_Mprime[row.names(mydf_Mprime) %in% dep_2$accession, ]
# row.names(norm_dep) <- norm_dep$Protein.Name
# heatmap_data <- as.data.frame(lapply(norm_dep[,5:22], as.numeric))
# heatmap_data <- as.data.frame(lapply(norm_dep[,c(5:10,13:22)], as.numeric))
# for (i in 1:length(row.names(heatmap_data))) {
#   row.names(heatmap_data) <- car::recode(row.names(heatmap_data), "dep_2$accession[i] = dep_2$Gene[i]")
# }

## without RRMS outlier
# heatmap_data <- heatmap_data[,c(1,2,5,6,11,12,15,16,23,26,7,8,13,14,19,20,3,4,9,10,17,18,21,22,24,25)] # with HC
# heatmap_data <- heatmap_data[,c(7,8,13,14,19,20,3,4,9,10,17,18,21,22,24,25)] # without HC
## no other outlier
# heatmap_data <- heatmap_data[,c(1,2,5,6,11,15,26,7,8,13,14,19,20,3,4,9,10,17,18,21,22,24,25)] # with HC
# heatmap_data <- heatmap_data[,c(6,11,12,15,16,23,26,7,8,13,14,19,20,3,4,10,17,18,22,24,25)] # with HC
# heatmap_data <- heatmap_data[,c(6,11,12,15,16,23,26,3,4,10,17,18,22,24,25)] # without RRMS

# progMS_HC <- read.xlsx("SPMS_vs_HC_stats_lowess_sex_ANCOVA_imputation.xlsx")
# progMS_HC_up <- progMS_HC[progMS_HC$diff_exp > 0 & progMS_HC$p_val < 0.1,]
# dep_2_up <- dep_2[dep_2$diff_exp > 0,]
# common.list_up <- Reduce(intersect, list(dep_2_up$accession,progMS_HC_up$accession))
# progMS_HC_down <- progMS_HC[progMS_HC$diff_exp < 0 & progMS_HC$p_val < 0.1,]
# dep_2_down <- dep_2[dep_2$diff_exp < 0,]
# common.list_down <- Reduce(intersect, list(dep_2_down$accession,progMS_HC_down$accession))
# common.list <- c(common.list_down,common.list_up)
# 
# targets <- t(read.delim("targets.txt", header = F))
# targets <- t(read.delim("ProgMS_HC_targets.txt", header = F))
# 
# heatmap_data <- mydf_Mprime[row.names(mydf_Mprime) %in% targets, ]

annotation_col <- data.frame(Lesion6 = sampleTable_sub$Lesion_type_6,
                             # Lesion9 = sampleTable_sub$Lesion_type_9,
                             sex = sampleTable_sub$sex,
                             row.names = sampleTable_sub$sampleName)
annotation_col$Lesion6 <- factor(annotation_col$Lesion6, levels = c("CWM","NAWM","2","3","4","6"))
annotation_col <- annotation_col[order(annotation_col$Lesion6),]
new_order <- rownames(annotation_col)
heatmap_data <- heatmap_data[, new_order]
# heatmap_data <- heatmap_data[targets,]
# heatmap_data <- heatmap_data[,11:26]
# row.names(heatmap_data) <- sub("(.*)_[0-9]$","\\1", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("_", "/", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("FA", "arachidonate", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("AcCa", "Acylcarnitine", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("LPC", "Lysophosphatidylcholine", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("PC", "Phosphatidylcholine", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("LPC", "Lysophosphatidylcholine", row.names(heatmap_data))
# row.names(heatmap_data) <- gsub("PC", "Phosphatidylcholine", row.names(heatmap_data))


# annotation_col <- annotation_col["disease"]

paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
heat <- pheatmap::pheatmap(heatmap_data, 
                           scale = "row", 
                           annotation_col = annotation_col,
                           cluster_cols = F, 
                           cluster_rows = T,
                           # clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                           color = colors, 
                           show_rownames = T, 
                           show_colnames = F,
                           angle_col = 45,
                           # fontsize_row = 5
                           #, cellwidth = 12   # decrease for more compressed plot
                           # ,cellheight = 10
)

## heatmap with target genes highlighted

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

targets <- c("PC(16:0_16:0)",
             "PC(16:0_20:3)",
             "PC(16:0_20:4)",
             "PC(16:0_22:4)",
             "PC(16:0_22:5)",
             "PC(16:0_22:6)",
             "PC(16:0, 20:4)_2",
             "PC(16:1_20:4)",
             "PC(18:0_20:4)",
             "PC(18:0_22:5)",
             "PC(18:0_22:6)",
             "PC(18:1_22:5)",
             "LPC(17:1)",
             "LPC(18:0)",
             "LPC(18:2)",
             "LPC(20:1)",
             "LPC(22:0)",
             "LPC(22:1)",
             "LPC(23:0)",
             "LPC(24:0)",
             "LPC(24:1)",
             "LPC(26:0)",
             "LPC(26:1)"
) # 

# hav2 <- grep("^PC\\(O-.*", dep_2$accession, value = T)

hav2 <- targets

add.flag(heat, kept.labels = hav2, repel.degree = 0.3)


### volcano

# sel1 <- my_res$Protein[my_res$stat == TRUE]
sel1 <- hav2

EnhancedVolcano(my_res,
                lab = my_res$accession, #NA,
                selectLab = sel1,#tmarkers[1,],
                x = 'diff_exp',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                xlim = c(-2, 3),
                ylim = c(0, 5),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 4.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "ProgMS vs. RRMS",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                #maxoverlapsConnectors = NULL,
                #min.segment.length = 0.0000001,
                #directionConnectors = "x",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)


## volcano 2 (under construction)

res <- my_res[!is.na(my_res$p_val),]
res$Gene <- res$accession
# selected_Proteins <- c("BDNF","TPPP3","CONA1")
selected_Proteins <- targets
# selected_Proteins <- my_res$Gene[my_res$stat == TRUE]

res$neglog10p <- -log10(res$p_val)

# Assign color categories
res$color <- "Not significant"
res$color[res$p_val < 0.05 & res$diff_exp > 0.1] <- "Upregulated"
res$color[res$p_val < 0.05 & res$diff_exp < -0.1] <- "Downregulated"

# Get coordinates for selected proteins
label_data <- res[res$Gene %in% selected_Proteins, ]

# Separate label data by upregulated and downregulated
label_data_up <- label_data[label_data$color == "Upregulated", ]
label_data_down <- label_data[label_data$color == "Downregulated", ]

# Assign x-labels for left (downregulated) and right (upregulated)
x_label_left <- min(res$diff_exp) - 1.5
x_label_right <- max(res$diff_exp) + 1.5

# Stacked y-labels for each side
y_up <- seq(max(res$neglog10p), min(res$neglog10p), length.out = max(1, nrow(label_data_up)))
y_down <- seq(max(res$neglog10p), min(res$neglog10p), length.out = max(1, nrow(label_data_down)))

# Data frames for label positions
stack_labels_up <- data.frame(
  Gene = label_data_up$Gene,
  x_label = x_label_right,
  y_label = y_up
)
stack_labels_down <- data.frame(
  Gene = label_data_down$Gene,
  x_label = x_label_left,
  y_label = y_down
)

# Merge to get original points
stack_labels_up <- merge(stack_labels_up, label_data_up, by = "Gene")
stack_labels_down <- merge(stack_labels_down, label_data_down, by = "Gene")
stack_labels <- rbind(stack_labels_up, stack_labels_down)

# Volcano plot
ggplot(res, aes(x = diff_exp, y = neglog10p)) +
  geom_point(aes(color = color), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Not significant" = "grey60", "Downregulated" = "lightblue")) +
  geom_segment(
    data = stack_labels,
    aes(x = diff_exp, y = neglog10p, xend = x_label, yend = y_label),
    color = "gray50", linewidth = 0.4, inherit.aes = FALSE
  ) +
  geom_label(
    data = stack_labels,
    aes(x = x_label, y = y_label, label = Gene),
    hjust = ifelse(stack_labels$x_label > 0, 0.65, 0.1), # right for upregulated, left for downregulated
    fontface = "bold", fill = "white", color = "black", label.size = 0.1, inherit.aes = FALSE
  ) +
  xlab(expression("Log"[2]*"(Fold-Change)")) +
  ylab(expression("-Log"[10]*"("*italic("p")*"-value)")) +
  ggtitle("MS Lesion vs. CWM - lipidomics") +
  coord_cartesian(xlim = c(x_label_left - .0, x_label_right + .6)) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())



#### GSEA
library(fgsea)
library(msigdbr)
library(dplyr)

deg_df <- my_res
deg_df$gene <- deg_df$accession
# colnames(deg_df)[1] <- "gene"
deg_df$log2FC <- deg_df$diff_exp
deg_df$pval <- deg_df$p_val

# Basic checks
stopifnot(all(c("gene", "log2FC", "pval") %in% colnames(deg_df)))

## 2. Build ranked gene list ----------------------------------------------
# Use signed -log10(p) so more significant genes with positive FC are high,
# and significant down-genes are strongly negative.[web:3][web:7][web:8]
deg_ranked <- deg_df %>%
  dplyr::filter(!is.na(pval), !is.na(log2FC)) %>%
  dplyr::mutate(rank_metric = sign(log2FC) * (-log10(pval))) %>%
  dplyr::arrange(desc(rank_metric))

# Create named numeric vector for fgsea
gene_list <- deg_ranked$rank_metric
names(gene_list) <- deg_ranked$gene

# Optional: remove duplicates, keep the most extreme rank per gene
gene_list <- tapply(gene_list, names(gene_list), function(x) x[which.max(abs(x))])
gene_list <- sort(unlist(gene_list), decreasing = TRUE)

## 3. Get gene sets (example: Hallmark from MSigDB) -----------------------
# msigdbr provides MSigDB gene sets directly in R.[web:3][web:7]
species <- "Homo sapiens"
msigdb_cat <- "C2"   # "H" = Hallmark; others: "C2", "C5", etc.

msig_tidy <- msigdbr(species = species, collection = msigdb_cat)

# Convert to list: pathway_name -> vector of genes
pathways_list <- msig_tidy %>%
  split(.$gs_name) %>%
  lapply(function(df) df$gene_symbol)

## 4. Run GSEA with fgsea -------------------------------------------------
set.seed(123)
fgsea_res <- fgsea(
  pathways = pathways_list,
  stats     = gene_list,
  minSize   = 15,      # tune as desired
  maxSize   = 500,     # tune as desired
  nperm     = 10000    # increase for more stable p-values
)

# Order by adjusted p-value
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

## 5. Save and inspect results --------------------------------------------
# Save table
write.xlsx(as.data.frame(fgsea_res), rowNames = T,file="GSEA_Reactome_RNAseq_Lesion6_IL_vs_CWM_stats_lowess_sex_age_ANCOVA_imputation_results.xlsx")

# Look at top enriched pathways
head(fgsea_res[, c("pathway", "NES", "pval", "padj")])

## 6. Example: plot enrichment curve for one pathway ----------------------
top_pathway <- fgsea_res$pathway[1]
plotEnrichment(pathways_list[[top_pathway]], gene_list) +
  ggplot2::ggtitle(top_pathway)


#####################################  pathways

file1 <- read.delim("pathways_RNAseq_Lesion6_IL_vs_CWM_stats_lowess_sex_age_ANCOVA_imputation_results_sel2.csv", header = T)
colunas1 <- file1$Term[1:23]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3.3,3.9,4.5)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=14),
    #axis.text.x = element_text(face="bold", color = "black",size=12), 
    plot.title = element_text(face="bold",color = "black",size=16, hjust = 1),
    axis.text.y = element_text(face="bold",color = "black",size=16),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Upregulated DEGs enriched pathways - Inactive lesion vs. CWM") +
  # theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()




### correlation 

### SDC2 accession is P34741

gene1 <- "FLVCR2"
gene1 <- "LPCAT1"
SDC2_exp <- t(mydf_Mprime[gene1,])
colnames(SDC2_exp) <- c(gene1)
SDC2_exp2 <- as.data.frame(SDC2_exp)

SDC2_exp2$sampleName <- row.names(SDC2_exp2)
SDC2_exp2 <- SDC2_exp2 %>%
  left_join(sampleTable_sub[, c("sampleName","sex","age","Lesion_type_9","Lesion_type_6")], by = "sampleName")


coluna1 <- c("CWM","NAWM","2","3","4","6")
my_comparisons <- list(c("CWM","NAWM"),
                       c("CWM","2"),
                       c("CWM","3"),
                       c("CWM","4"),
                       c("CWM","6"),
                       c("NAWM","2"),
                       c("NAWM","3"), 
                       c("NAWM","4"), 
                       c("NAWM","6"),
                       c("2","3"),
                       c("2","4"),
                       c("2","6"),
                       c("3","4"),
                       c("3","6"),
                       c("4","6")
                       )

coluna1 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6")

p1 <- ggplot(SDC2_exp2, aes_string(x="factor(Lesion_type_6,levels=coluna1)", y=gene1, color="Lesion_type_6")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  # scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("Normalized gene expression") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle(gene1) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(face="bold", size=14),
                                                        axis.text.x = element_text(face="bold", size=16, colour = "black"),
                                                        axis.text.y = element_text(face="bold", size=14, colour = "black"))
p1











### replicate average
SDC2_exp3 <- data.table(SDC2_exp2, key = c("barcode"))
SDC2_exp3 <- SDC2_exp3[, list(SDC2=mean(SDC2)), by=c("barcode","Treatment","LME","CL.vol","CL.num")]
SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$Treatment %in% c("post")),]
# SDC2_exp3 <- SDC2_exp3[SDC2_exp3$SDC2 > -3,]
SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$LME %in% c("7")),]
# SDC2_exp3 <- SDC2_exp3[c(1,2,4:8),]

SDC2_exp3 <- SDC2_exp2[c(1:4,6:16,18),]
SDC2_exp3 <- SDC2_exp2[c(9:16,18),]
SDC2_exp3 <- SDC2_exp2[c(1:4,6:14,16),]
SDC2_exp3 <- SDC2_exp2

SDC2_exp3 <- SDC2_exp2[SDC2_exp2$SDC2 > -6,]
SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$Treatment %in% c("post")),]
SDC2_exp3 <- SDC2_exp3[c(3:6,9:38),]

# ggscatter(SDC2_exp3, x = "SDC2", y = "walktime", #label = rownames(mydata3),
#           repel = T, #font.label = c(14, "bold.italic", "black"),
#           add = "reg.line", # Add regression line
#           add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
#           conf.int = T ) + # Add confidence interval
#   stat_cor(method = "pearson", label.x = 0.5, label.y = 20) +
#   annotate("text", x = 3, y = 22, label = (paste0("slope==", coef(lm(SDC2_exp3$walktime~SDC2_exp3$SDC2))[2])), parse = T) + 
#   ggtitle(expression("")) + ylab(expression("25 ft walk time")) +
#   xlab(expression(italic("SDC2")*" normalized expression"))

# ggscatter(SDC2_exp3, x = "SDC2", y = "EDSS", #label = rownames(mydata3),
#           repel = T, #font.label = c(14, "bold.italic", "black"),
#           add = "reg.line", # Add regression line
#           add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
#           conf.int = T ) + # Add confidence interval
#   stat_cor(method = "pearson", label.x = 0.5, label.y = 20) +
#   annotate("text", x = 3, y = 22, label = (paste0("slope==", coef(lm(SDC2_exp3$EDSS~SDC2_exp3$SDC2))[2])), parse = T) + 
#   ggtitle(expression("")) + ylab(expression("EDSS")) +
#   xlab(expression(italic("SDC2")*" normalized expression"))

ggscatter(SDC2_exp3, x = "SDC2", y = "CL.num", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = -0.8, label.y = 27) +
  annotate("text", x = -0.6, y = 30, label = (paste0("slope==", coef(lm(SDC2_exp3$CL.num~SDC2_exp3$SDC2))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("CL number")) +
  xlab(expression(italic("SDC2")*" normalized expression"))

ggscatter(SDC2_exp3, x = "SDC2", y = "CL.vol", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = -0.8, label.y = .7) +
  annotate("text", x = -0.6, y = .6, label = (paste0("slope==", coef(lm(SDC2_exp3$CL.vol~SDC2_exp3$SDC2))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("CL Volume (mL)")) +
  xlab(expression(italic("SDC2")*" normalized expression"))

ggscatter(SDC2_exp3, x = "SDC2", y = "LME", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = -3.2, label.y = 6.5) +
  annotate("text", x = -2.5, y = 6, label = (paste0("slope==", coef(lm(SDC2_exp3$LME~SDC2_exp3$SDC2))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("LME foci")) +
  xlab(expression("CD362 smoothed expression (LOWESS)"))

# ggscatter(SDC2_exp3, x = "SDC2", y = "LME", #label = rownames(mydata3),
#           repel = T, #font.label = c(14, "bold.italic", "black"),
#           add = "reg.line", # Add regression line
#           add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
#           conf.int = T ) + # Add confidence interval
#   stat_cor(method = "spearman", label.x = -2., label.y = 7.8) +
#   annotate("text", x = -1.5, y = 7, label = (paste0("slope==", coef(lm(SDC2_exp3$LME~SDC2_exp3$SDC2))[2])), parse = T) + 
#   ggtitle(expression("RRMS and ProgMS")) + ylab(expression("LME foci")) +
#   xlab(expression("CD362 smoothed expression (LOWESS)"))
# 
# 
# 
# ### LME only ProgMS
# ggscatter(SDC2_exp3, x = "SDC2", y = "LME", #label = rownames(mydata3),
#           repel = T, #font.label = c(14, "bold.italic", "black"),
#           add = "reg.line", # Add regression line
#           add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
#           conf.int = T ) + # Add confidence interval
#   stat_cor(method = "pearson", label.x = -0.5, label.y = 8.8) +
#   annotate("text", x = -0.27, y = 8, label = (paste0("slope==", coef(lm(SDC2_exp3$LME~SDC2_exp3$SDC2))[2])), parse = T) + 
#   ggtitle(expression("ProgMS samples")) + ylab(expression("LME foci")) +
#   xlab(expression("CD362 smoothed expression (LOWESS)"))
# 
# ## LME all samples
# ggscatter(SDC2_exp3, x = "SDC2", y = "LME", #label = rownames(mydata3),
#           repel = T, #font.label = c(14, "bold.italic", "black"),
#           add = "reg.line", # Add regression line
#           add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
#           conf.int = T ) + # Add confidence interval
#   stat_cor(method = "pearson", label.x = -0.8, label.y = 7.8) +
#   annotate("text", x = -0.5, y = 7, label = (paste0("slope==", coef(lm(SDC2_exp3$LME~SDC2_exp3$SDC2))[2])), parse = T) + 
#   ggtitle(expression("RRMS and ProgMS")) + ylab(expression("LME foci")) +
#   xlab(expression("CD362 smoothed expression (LOWESS)"))
# 



### replicate average
SDC2_exp3 <- data.table(SDC2_exp2, key = c("barcode"))
SDC2_exp3 <- SDC2_exp3[, list(SDC2=mean(SDC2)), by=c("barcode","Treatment","LME","CL.vol","CL.num")]
SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$Treatment %in% c("post")),]

SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$Treatment %in% c("post")),]
SDC2_exp3 <- SDC2_exp3[SDC2_exp3$SDC2 > -3,]
SDC2_exp3 <- SDC2_exp3[!(SDC2_exp3$LME %in% c("7")),]
SDC2_exp3$LMEcat <- with(SDC2_exp3, ifelse(LME > 2, "High", "Low"))

coluna1 <- c("Low","High")
my_comparisons <- list(c("Low","High"))

p1 <- ggplot(SDC2_exp3, aes_string(x="factor(LMEcat,levels=coluna1)", y="SDC2", color="LMEcat")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  # ggtitle("COX2+ Classical monocytes") +
  # ggtitle("Nonclassical monocytes") +
  ggtitle("") +
  # ggtitle("Intermediate monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1


SDC2_exp3$LMEcat <- factor(SDC2_exp3$LMEcat, levels = c("Low","High"))
colors <- brewer.pal(n=4,name = "Dark2")
colors <- colors[c(1,3)]
plt <- ggbetweenstats(data = SDC2_exp3, x = LMEcat, y = SDC2,
                      # plot.type = "violin",
                      p.adjust.method = "none", 
                      xlab = "", ylab = "",
                      ggplot.component = list(theme(axis.text.x = element_text(size = 14, face = "bold"),
                                                    axis.text.y = element_text(size = 14, face = "bold"),
                                                    axis.title.y = element_text(size = 14, face = "bold"))),
                      ggsignif.args = list(textsize = 4),
                      boxplot.args = list(width = 0),
                      point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 1, size = 3, stroke = 0),
                      violin.args = list(width = 0.5, alpha = 0.2),
                      package = "RColorBrewer", palette = "Dark2",
                      results.subtitle = FALSE,
                      mean.plotting = FALSE,
                      mean.point.args = list(size = 0, color = "darkred"),
                      mean.label.args = list(size = 0),
                      messages = FALSE,
                      # package = "ggsci", palette = "nrc_npg", 
                      # package = "ggsci", palette = "uniform_startrek",
                      # package = "yarrr", palette = "info2",
                      type = "np"
)
plt


















