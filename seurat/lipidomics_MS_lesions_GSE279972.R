#############################################
####################
########## starting from raw abundance

############## lipidomics


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



setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/bulkRNAseq_MS_lesions_GSE279972/")
# setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/MS/lipidomics_monocytes/no_normalization/")

data1 <- read.xlsx("Processed_data_all_omics.xlsx", sheet = 3)
data1$Accession <- data1$X1
data1_unique <- data1[!duplicated(data1$Accession), ]

dat0 <- data1_unique[,c(2:101)] # all samples
# dat <- dat[,c(3:25,28:30)] # outlier removal
dat <- 2^dat0
row.names(dat) <- data1_unique$Accession
dat_numeric <- as.data.frame(lapply(dat, as.numeric))
row.names(dat_numeric) <- data1_unique$Accession

stopifnot(is.matrix(dat_numeric) || is.data.frame(dat_numeric))
dat_numeric <- as.matrix(dat_numeric)
## keep features with >=60% non-missing overall
## and in each biological group (here only overall, an extend by group)
na_prop <- rowMeans(is.na(dat_numeric))
keep_features <- na_prop <= 0.40        # <=40% missing
dat_numeric_filt <- dat_numeric[keep_features, , drop = FALSE]
## Random forest imputation for mostly MCAR/MAR missingness
set.seed(1)
rf_res <- missForest(t(dat_numeric_filt))   # missForest expects samples in rows
dat_numeric_rf_imp <- t(rf_res$ximp)        # back to features in rows
# ### QRILC is recommended for left-censored MNAR; log-transform first 
# log_dat_numeric <- log2(dat_numeric_filt)          # or log1p() if zeros are present
# qr_res <- impute.QRILC(log_dat_numeric)    # returns a list: [[1]] = imputed matrix 
# log_dat_numeric_qrilc <- qr_res[[1]]
# dat_numeric_qrilc_imp <- 2^log_dat_numeric_qrilc   # back-transform to original scale


# dat_numeric <- dat_rf_imp
# dat_numeric <- dat
# dat_numeric <- as.data.frame(lapply(dat, as.numeric))
# row.names(dat_numeric) <- data1_unique$Accession
#### more stringent
raw0 <- dat_numeric_rf_imp[complete.cases(dat_numeric_rf_imp) & apply(dat_numeric_rf_imp != 0, 1, all),]
raw <- as.data.frame(raw0)

## norm log2(value/colMedian)
mydfmedian <- raw / colMedians(data.matrix(raw))[col(raw)]
lg2mydfmean <- log2(mydfmedian)

mydf_M <- lg2mydfmean - lg2mydfmean[,1]
mydf_A <- lg2mydfmean
# mydf_A <- data.frame(matrix(ncol = 18,nrow = 2262))

for (i in 1:ncol(lg2mydfmean)) {
  # avg_df[[paste0("avg_", names(lg2mydfmean)[i])]] <- rowMeans(cbind(lg2mydfmean[[1]], lg2mydfmean[[i]]), na.rm = TRUE)
  mydf_A[[names(lg2mydfmean)[i]]] <- rowMeans(cbind(lg2mydfmean[[1]], lg2mydfmean[[i]]), na.rm = TRUE)
}

window_size <- 40

# Initialize output dataframe to store moving averages for M columns
mydf_M_avg <- data.frame(matrix(NA, nrow = nrow(mydf_M), ncol = ncol(mydf_M)))
colnames(mydf_M_avg) <- paste0("M_avg_", colnames(mydf_M))

for (i in 1:ncol(mydf_A)) {
  # Combine A and M columns
  temp_df <- data.frame(A = mydf_A[[i]], M = mydf_M[[i]])
  # Sort by A ascending
  temp_df_sorted <- temp_df[order(temp_df$A), ]
  # Calculate moving average on M (sorted by A)
  temp_df_sorted$M_avg <- rollmean(temp_df_sorted$M, k = window_size, fill = NA, align = "right")
  # Replace NA in M_avg with original M values at those positions
  nas <- is.na(temp_df_sorted$M_avg)
  temp_df_sorted$M_avg[nas] <- temp_df_sorted$M[nas]
  # Reorder M_avg back to original row order by using order index
  # Create a vector to map sorted indices to original locations
  original_order <- order(order(temp_df$A))
  # Put the moving average back to the original dataframe row order
  mydf_M_avg[[i]] <- temp_df_sorted$M_avg[original_order]
}


row.names(mydf_M_avg) <- row.names(mydf_M)
mydf_Mprime <- mydf_M - mydf_M_avg


meta <- read.xlsx("Processed_data_all_omics.xlsx", sheet = 2)
meta$ID0 <- gsub(" ",".", meta$Unifying_code)
meta$ID <- gsub("-",".", meta$ID0)
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

### t test
######## run t test
# ## Logical indices for "RRMS","SPMS","HC" samples
# idx_prog  <- sampleTable_sub$group == "ProgMS"
# idx_st    <- sampleTable_sub$group == "RRMS"
# 
# ## Container for results
# res_ttest <- data.frame(
#   metabolite = rownames(mydf_Mprime),
#   diff_exp = NA_real_,
#   pvalue_ttest = NA_real_,
#   stringsAsFactors = FALSE
# )
# 
# ## Loop over metabolites, compute mean difference and t-test p-value
# for (i in seq_len(nrow(mydf_Mprime))) {
#   y <- as.numeric(mydf_Mprime[i, ])
#   ## Expression vectors per group
#   y_prog <- y[idx_prog]
#   y_st   <- y[idx_st]
#   ## Remove NAs
#   y_prog <- y_prog[!is.na(y_prog)]
#   y_st   <- y_st[!is.na(y_st)]
#   ## Need at least 2 per group for t-test to be meaningful
#   if (length(y_prog) < 2 || length(y_st) < 2) next
#   ## Mean difference ProgMS - StRRMS
#   res_ttest$diff_exp[i] <- mean(y_prog) - mean(y_st)
#   ## Two-sample t-test (by default, Welch's t-test)
#   tt <- t.test(y_prog, y_st)
#   res_ttest$pvalue_ttest[i] <- tt$p.value
# }
# 
# my_res <- my_res[order(my_res$p_val),]
# my_res$index <- 1:nrow(my_res)
# alpha <- 0.5
# my_res_length <- length(my_res$index)
# my_res$FDR <- my_res$index*alpha / my_res_length
# my_res$stat <- with(my_res, ifelse(p_val < FDR, TRUE, FALSE))

########

my_res <- res

# my_res$Gene <- my_res$accession
# my_res$Protein <- my_res$accession
# 
# for (i in 1:length(my_res$Protein)) {
#   my_res$Gene <- car::recode(my_res$Gene, "data1_unique$Accession[i] = data1_unique$Genes[i]")
#   # my_res$Protein <- car::recode(my_res$Protein, "data1_unique$Accession[i] = data1_unique$Protein.Name[i]")
# }

# Export results
write.xlsx(as.data.frame(my_res), rowNames = T, file="lipidomics/Lesion6_IL_vs_CWM_stats_lowess_sex_age_ANCOVA_imputation.xlsx")




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


### correlation 

### SDC2 accession is P34741

SDC2_exp <- t(mydf_Mprime["P34741",])
colnames(SDC2_exp) <- c("SDC2")
# colnames(SDC2_exp) <- sub("(.*)_.*_.*_.*","\\1",colnames(dds3_CM))
# rownames(SDC2_exp)[1] <- "SDC2"
# SDC2_exp2 <- t(SDC2_exp)
SDC2_exp2 <- as.data.frame(SDC2_exp)
SDC2_exp2$CL.vol <- row.names(SDC2_exp2)
SDC2_exp2$CL.num <- row.names(SDC2_exp2)
SDC2_exp2$LME <- row.names(SDC2_exp2)
SDC2_exp2$Treatment <- row.names(SDC2_exp2)
# SDC2_exp2$EDSS <- row.names(SDC2_exp2)
# SDC2_exp2$walktime <- row.names(SDC2_exp2)
SDC2_exp2$barcode <- row.names(SDC2_exp2)

for (i in 1:length(sampleTable$sampleName)){
  SDC2_exp2$CL.num <- car::recode(SDC2_exp2$CL.num, "sampleTable$sampleName[i] = sampleTable$CL.num[i]")
  SDC2_exp2$CL.vol <- car::recode(SDC2_exp2$CL.vol, "sampleTable$sampleName[i] = sampleTable$CL.vol[i]")
  SDC2_exp2$LME <- car::recode(SDC2_exp2$LME, "sampleTable$sampleName[i] = sampleTable$LME[i]")
  # SDC2_exp2$EDSS <- car::recode(SDC2_exp2$EDSS, "meta$GEM[i] = meta$EDSS[i]")
  # SDC2_exp2$walktime <- car::recode(SDC2_exp2$walktime, "meta$GEM[i] = meta$X25.FT.Walk.Time..s.[i]")
  SDC2_exp2$barcode <- car::recode(SDC2_exp2$barcode, "sampleTable$sampleName[i] = sampleTable$patient[i]")
  SDC2_exp2$Treatment <- car::recode(SDC2_exp2$Treatment, "sampleTable$sampleName[i] = sampleTable$treatment[i]")
}

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


















