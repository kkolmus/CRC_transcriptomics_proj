### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(GEOquery)
  library(tidyverse)
})

### MAIN ###

flog.threshold(NULL)

flog.debug("Prepare working environment")

rm(list = ls())

proj.dir <- "~/Desktop/HT projects/CRC_transcriptomics_proj"
data.dir <- file.path(proj.dir, "data")
input.dir <- file.path(data.dir, "input")
dir.create(input.dir, recursive = TRUE)
output.dir <- file.path(data.dir, "output")
dir.create(output.dir, recursive = TRUE)

flog.debug("Download metadata for GSE50760")

projID <- "GSE50760"
temp.dir <- file.path(input.dir, projID)
dir.create(temp.dir, recursive = TRUE)

gprojID <- getGEO(
  GEO = projID,
  destdir = temp.dir,
  GSEMatrix = TRUE)


flog.debug("Clean metadata")

pdata <- read.table(file.path(input.dir, projID, "GSE50760_series_matrix.txt"), fill = TRUE, skip = 31)
pdata <- as.data.frame(t(pdata))
temp_colnames <- unname(unlist(pdata[1,, drop <- TRUE]))
temp_colnames <- as.character(temp_colnames)
temp_colnames <- substring(temp_colnames, 9, nchar(temp_colnames))
colnames(pdata) <- temp_colnames
pdata <- pdata[-1,-c(39:41)]
pdata <- pdata[, c(2,1,3:ncol(pdata))]
rownames(pdata) <- NULL
rownames(pdata) <- pdata$geo_accession

pdata <- pdata[, c(1,2,8,10,11)]
new_colnames <- c("geo_accession", "patient_id", "source_tissue", "diagnosis", "ajcc_stage")
colnames(pdata) <- new_colnames
pdata$patient_id <- gsub(".*\\s", "", pdata$patient_id)
pdata$diagnosis <- substring(pdata$diagnosis, 9, 100)
pdata$ajcc_stage <- substring(pdata$ajcc_stage, 13, 100)

pdata <- filter(pdata, pdata$source_tissue != "metastasized cancer")

saveRprojID(pdata, file.path(temp.dir, "pdata.RprojID"))


flog.debug("Download normalized count data for GSE50760")

getGEOSuppFiles(
  GEO = projID,
  baseDir = temp.dir
)

setwd(file.path(temp.dir, projID))
getwd()

files <- list.files(path = file.path(temp.dir, projID))
untar(tarfile <- file.path(files))
file.remove(file.path(files))

files <- list.files(path = file.path(temp.dir, projID))
for (f in files) {
  gunzip(f)
}

files <- list.files(path = file.path(temp.dir, projID))
for (f in files) {
  temp_file_name = substring(f, 1, 10)
  temp_file = read.table(file.path(temp.dir, projID, f), header = TRUE)
  temp_file$genes = as.character(temp_file$genes)
  temp_file$genes = make.unique(temp_file$genes)
  temp_file = column_to_rownames(temp_file, "genes")
  assign(temp_file_name, temp_file)
}

pattern <- grep("GSM", names(.GlobalEnv), value = TRUE)
listDF <- do.call("list",mget(pattern))
rm(list = ls(pattern = "GSM"))

n.counts <- listDF[[1]]
for (i in 2:length(listDF)) {
  n.counts = merge(n.counts, listDF[[i]], by = "row.names", all = TRUE)
  rownames(n.counts) = n.counts$Row.names
  n.counts = n.counts[ , !(names(n.counts) %in% "Row.names")]  
}


flog.debug("Clean normalized count data for GSE50760")

temp_colnames <- colnames(n.counts)
temp_colnames <- gsub("_FPKM", "", temp_colnames)
temp_colnames <- gsub("\\.", "-", temp_colnames)

colnames(n.counts) <- temp_colnames

names2use <- pdata$patient_id
n.counts <- n.counts[, names(n.counts) %in% names2use]

stopifnot(pdata$patient_id %in% colnames(n.counts))

n.counts <- rownames_to_column(n.counts, "Symbol")

saveRDS(n.counts, file.path(temp.dir, "n.counts.RDS"))