### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(GEOquery)
  library(ArrayExpress)
  library(tidyverse)
})

`%!in%` <- Negate(`%in%`)

### MAIN ###

flog.threshold(NULL)

flog.debug("Prepare working environment")

rm(list = ls())

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
data.dir = file.path(proj.dir, "data")
input.dir = file.path(data.dir, "input")
dir.create(input.dir, recursive = TRUE)
output.dir = file.path(data.dir, "output")
dir.create(output.dir, recursive = TRUE)

# load RNA-Seq and microarray dataset
source(file = file.path(data.dir, "projects.R"))
RNA_Seq_datasets <- RNA_Seq_datasets[-5]

flog.debug("Download datasets")

for (ds in RNA_Seq_datasets) {
  flog.debug(paste0("Downloading data for the ", ds, " dataset"))
  
  temp.dir = file.path(input.dir, ds)
  dir.create(temp.dir, recursive = TRUE)
  
  if (grepl(pattern = "GSE", x = ds) == TRUE) {
    gds = getGEO(
      GEO = ds,
      destdir = temp.dir,
      GSEMatrix = TRUE)
    getGEOSuppFiles(
      GEO = ds,
      baseDir = temp.dir
      )
    # pheno data
    pdata = gds[[1]]@phenoData@data
    assign(paste0("pdata_", ds), pdata)
  } else {
     anno_AE = getAE(
       accession = ds,
       path = temp.dir, 
       type = "processed")
    }
}


flog.debug("Analyze metadata")

###############################
flog.debug("Analyze GSE137327")
###############################

pdata_GSE137327 <- pdata_GSE137327[, c(1,2,39,40)] 
pdata_GSE137327 <- pdata_GSE137327[, c(2,1,4,3)]
new_colnames <- c("geo_accession", "patient_id", "source_tissue", "diagnosis")
colnames(pdata_GSE137327) <- new_colnames
pdata_GSE137327$status <- rep(c("Cancer", "Colon"), 9)

saveRDS(pdata_GSE137327, file.path(input.dir, "GSE137327", "pdata.RDS"))

setwd(file.path(input.dir, "GSE137327", "GSE137327"))
getwd()
untar(tarfile = file.path(input.dir, "GSE137327", "GSE137327", "GSE137327_RAW.tar"))
file.remove(file.path(input.dir, "GSE137327", "GSE137327", "GSE137327_RAW.tar"))

files <- list.files(path = file.path(input.dir, "GSE137327", "GSE137327"))
for (f in files) {
  gunzip(f)
}

files <- list.files(path = file.path(input.dir, "GSE137327", "GSE137327"))
for (f in files) {
  temp_file_name = substring(f, 1, 10)
  temp_file = read.table(file.path(input.dir, "GSE137327", "GSE137327", f), fill = TRUE)
  temp_file = filter(temp_file, str_detect(temp_file$V2, "NM_"))
  temp_file = temp_file[, c(6,5)]
  colnames(temp_file) <- c("genes", paste0("FPKM_", temp_file_name))
  temp_file$genes = as.character(temp_file$genes)
  temp_file[,2] = as.numeric(as.character(temp_file[,2]))
  temp_file = na.omit(temp_file)
  temp_file$genes = make.unique(temp_file$genes)
  rownames(temp_file) <- NULL
  temp_file = column_to_rownames(temp_file, "genes")
  assign(temp_file_name, temp_file)
}

pattern <- grep("GSM", names(.GlobalEnv), value = TRUE)
listDF <- do.call("list", mget(pattern))
rm(list = ls(pattern = "GSM"))

n.counts <- listDF[[1]]
for (i in 2:length(listDF)) {
  n.counts = merge(n.counts, listDF[[i]], by = "row.names", all = TRUE)
  rownames(n.counts) = n.counts$Row.names
  n.counts = n.counts[ , !(names(n.counts) %in% "Row.names")]  
}

temp_colnames <- colnames(n.counts)
temp_colnames <- gsub("FPKM_", "", temp_colnames)
temp_colnames <- paste(temp_colnames, pdata_GSE137327$status, sep = "_")

colnames(n.counts) <- temp_colnames

n.counts <- rownames_to_column(n.counts, "Symbol")

saveRDS(n.counts, file.path(input.dir, "GSE137327", "n.counts.RDS"))


###############################
flog.debug("Analyze GSE144259")
###############################

pdata_GSE144259 <- pdata_GSE144259[,c(2,1,8)]
pdata_GSE144259$source_tissue <- rep(c("Colon", "Cancer", "Liver_metastasis"), 3)
pdata_GSE144259 <- filter(pdata_GSE144259, pdata_GSE144259$source_tissue %!in% "Liver_metastasis")
pdata_GSE144259 <- pdata_GSE144259[,c(1,2,4,3)]
new_colnames <- c("geo_accession", "patient_id", "source_tissue", "diagnosis")
colnames(pdata_GSE144259) <- new_colnames

saveRDS(pdata_GSE144259, file.path(input.dir, "GSE144259", "pdata.RDS"))

setwd(file.path(input.dir, "GSE144259", "GSE144259"))
getwd()

files <- list.files(path = file.path(input.dir, "GSE144259", "GSE144259"))
for (f in files) {
  gunzip(f)
}

files <- list.files(path = file.path(input.dir, "GSE144259", "GSE144259"))

n.counts <- read.table(file.path(files), header = TRUE)
n.counts$GeneID <- as.character(n.counts$GeneID)
n.counts$GeneID <- make.unique(n.counts$GeneID)
n.counts <- column_to_rownames(n.counts, "GeneID")

names2use <- pdata_GSE144259$patient_id
n.counts <- n.counts[, names(n.counts) %in% names2use]

n.counts <- rownames_to_column(n.counts, "Symbol")
n.counts <- filter_all(n.counts, all_vars(. != 0))

saveRDS(n.counts, file.path(input.dir, "GSE144259", "n.counts.RDS"))


###############################
flog.debug("Analyze GSE87096")
###############################

pdata_GSE87096 <- filter(pdata_GSE87096, pdata_GSE87096$molecule_ch1 == "total RNA")
pdata_GSE87096 <-  pdata_GSE87096[,c(2,1,39)]
pdata_GSE87096$title <- sub("(\\d)[^0-9]+$", "\\1", pdata_GSE87096$title)
pdata_GSE87096$source_tissue <- rep(c("Colon", "Cancer"), 6)
pdata_GSE87096 <-  pdata_GSE87096[,c(1,2,4,3)]

new_colnames <- c("geo_accession", "patient_id", "source_tissue", "diagnosis")

colnames(pdata_GSE87096) <- new_colnames

saveRDS(pdata_GSE87096, file.path(input.dir, "GSE87096", "pdata.RDS"))

setwd(file.path(input.dir, "GSE87096", "GSE87096"))
getwd()

file.remove("GSE87096_experiment.xml.gz", "GSE87096_RAW.tar", "GSE87096_run.xml.gz", "GSE87096_submission.xml.gz")

files <- list.files(path = file.path(input.dir, "GSE87096", "GSE87096"))
for (f in files) {
  gunzip(f)
}

files <- list.files(path = file.path(input.dir, "GSE87096", "GSE87096"))

for (f in files) {
  temp_file_name = substring(f, 1, 17)
  temp_file = read.table(file.path(input.dir, "GSE87096", "GSE87096", f), fill = TRUE, header = TRUE)
  temp_file = temp_file[, c(2,7,8)]
  temp_file = na.omit(temp_file)
  rownames(temp_file) = NULL
  temp_file = column_to_rownames(temp_file, "Gene_Name")
  assign(temp_file_name, temp_file)
}

pattern <- grep("GSE87096", names(.GlobalEnv), value = TRUE)
listDF <- do.call("list", mget(pattern))
rm(list = ls(pattern = "GSE87"))

n.counts <- listDF[[1]]
for (i in 2:length(listDF)) {
  n.counts = merge(n.counts, listDF[[i]], by = "row.names", all = TRUE)
  rownames(n.counts) = n.counts$Row.names
  n.counts = n.counts[ , !(names(n.counts) %in% "Row.names")]  
}

temp_colnames <- colnames(n.counts)
temp_colnames <- gsub(".RPKM", "", temp_colnames)
colnames(n.counts) <- temp_colnames

temp_colnames <- colnames(n.counts)
temp_colnames <- gsub("T", "_Cancer", temp_colnames)
temp_colnames <- gsub("N", "_Colon", temp_colnames)
colnames(n.counts) <- temp_colnames

n.counts <- rownames_to_column(n.counts, "Symbol")

saveRDS(n.counts, file.path(input.dir, "GSE87096", "n.counts.RDS"))


###############################
flog.debug("Analyze GSE92914")
###############################

pdata_GSE92914 <- pdata_GSE92914[, c(2,1,49,48)]
new_colnames <- c("geo_accession", "patient_id", "source_tissue", "diagnosis")
colnames(pdata_GSE92914) <- new_colnames

pdata_GSE92914$patient_id <- as.character(pdata_GSE92914$patient_id)

pdata_GSE92914$sampleID <- substring(
  pdata_GSE92914$patient_id, 
  (nchar(pdata_GSE92914$patient_id)-3), 
  nchar(pdata_GSE92914$patient_id))
pdata_GSE92914$sampleID <- paste0("X", pdata_GSE92914$sampleID)

pdata_GSE92914 <- pdata_GSE92914[-c(1,4,7),]
pdata_GSE92914 <- pdata_GSE92914[c(1:6),]

saveRDS(pdata_GSE92914, file.path(input.dir, "GSE92914", "pdata.RDS"))

setwd(file.path(input.dir, "GSE92914", "GSE92914"))
getwd()

files <- list.files(path = file.path(input.dir, "GSE92914", "GSE92914"))

untar("GSE92914_RAW.tar")

file.remove("GSE92914_RAW.tar")

files <- list.files(path = file.path(input.dir, "GSE92914", "GSE92914"))
for (f in files) {
  gunzip(f)
}

files <- list.files(path = file.path(input.dir, "GSE92914", "GSE92914"))

for (f in files) {
  temp_file_name = substring(f, 1, 17)
  temp_file = read.table(file.path(input.dir, "GSE92914", "GSE92914", f), fill = TRUE, header = TRUE)
  temp_file = na.omit(temp_file)
  rownames(temp_file) = NULL
  temp_file = column_to_rownames(temp_file, "Ensemble_ID")
  assign(temp_file_name, temp_file)
}

pattern <- grep("GSM", names(.GlobalEnv), value = TRUE)
listDF <- do.call("list", mget(pattern))
rm(list = ls(pattern = "GSM"))

n.counts <- listDF[[1]]
for (i in 2:length(listDF)) {
  n.counts = merge(n.counts, listDF[[i]], by = "row.names", all = TRUE)
  rownames(n.counts) = n.counts$Row.names
  n.counts = n.counts[ , !(names(n.counts) %in% "Row.names")]  
}

temp_colnames <- colnames(n.counts)
temp_colnames <- gsub("_RPKM", "", temp_colnames)
colnames(n.counts) <- temp_colnames

names2use <- pdata_GSE92914$sampleID
n.counts <- n.counts[, names(n.counts) %in% names2use]

n.counts <- rownames_to_column(n.counts, "ENSEMBLid")

saveRDS(n.counts, file.path(input.dir, "GSE87096", "n.counts.RDS"))

# concern regarding annotation of this dataset