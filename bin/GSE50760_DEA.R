### LIBS ###

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("apeglm")

suppressPackageStartupMessages({
  library(futile.logger)
  library(SummarizedExperiment)
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
})


### MAIN ###

flog.threshold(NULL)

flog.debug("Prepare working environment")

rm(list = ls())

projID =  "GSE50760"

proj.dir <- "~/Desktop/HT projects/CRC_transcriptomics_proj"
data.dir <- file.path(proj.dir, "data")
input.dir <- file.path(data.dir, "input", projID)
output.dir <- file.path(data.dir, "output", projID)


counts <- readRDS(file.path(input.dir, "counts.RDS"))
rownames(counts) <- NULL
counts <- column_to_rownames(counts, "ID")
counts <- counts[, -1]
counts <- as.matrix(counts)
pdata <- readRDS(file.path(input.dir, "pdata.RDS"))
pdata <- pdata[match(rownames(pdata), colnames(counts)),]
pdata$tissue_type <- ifelse(
  pdata$source_tissue == "primary colorectal cancer",
  yes = "cancer",
  no = "colon")

stopifnot(rownames(pdata) == colnames(counts))


flog.debug("Differential expression analysis")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = pdata,
  design = ~ tissue_type)
dds

# prefilter
reads = 50
# keep <- rowSums((counts(dds)) >= 25) >= (ncol((counts(dds))/2))
keep <- rowSums(counts(dds)) >= (ncol(counts(dds))*reads)
dds <- dds[keep,]
dds

# relevel
dds$tissue_type <- factor(dds$tissue_type, levels = c("colon","cancer"))
dds$tissue_type <- relevel(dds$tissue_type, ref = "colon")

# DEA
dds <- DESeq(dds)

res <- results(dds)
res

resultsNames(dds)


reference <- readRDS(file.path(input.dir, "counts.RDS"))[, c(1,2)]


res <- results(dds, name = "tissue_type_cancer_vs_colon")
res_gene <- res@rownames
res <- as.data.frame(res@listData)
rownames(res) <- res_gene
res <- rownames_to_column(res, "ID")

res <- full_join(res, reference, by = "ID")
res <- res[,c(1,8,2:7)]


raw_data <- as.data.frame(counts)
new_colnames <- paste0(pdata$geo_accession, "_", pdata$tissue_type)
colnames(raw_data) <- new_colnames
raw_data <- rownames_to_column(raw_data, "ID")


res <- left_join(res, raw_data, by = c("ID"))
res <- res[,-c(3,5:7)]


saveRDS(res, file.path(input.dir, "DEA.RDS"))