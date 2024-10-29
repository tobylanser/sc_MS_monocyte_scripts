setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper/compass")

library(data.table)

counts <- fread("compass_input/myeloid_counts.csv", nThread = 18)
gc()
counts <- counts[, -1]
gc()
counts[1,]

fwrite(counts, "compass_input/myeloid_counts.csv", sep = ",", row.names = T, col.names = T, nThread = 18)
