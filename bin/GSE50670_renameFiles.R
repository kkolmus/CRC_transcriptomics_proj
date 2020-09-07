### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
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
count.input.dir <- file.path(input.dir, paste0(projID, "_counts"))

files.ls <- list.files(path = count.input.dir, pattern = "count", full.names = FALSE)
samples <- substring(files.ls, 1, (nchar(files.ls)-6))

map <- read.csv(file.path(input.dir, "SraRunInfo.csv"))
map$Run <- as.character(map$Run)
map$SampleName <- as.character(map$SampleName)
map <- dplyr::filter(map, map$Run %in% samples)
map <- map[, c(1,30)]

for (i in 1:length(files.ls)) {
  temp_file = files.ls[[i]]
  file.rename(
    from = file.path(count.input.dir, temp_file), 
    to = file.path(count.input.dir, paste0(map[i,2], ".csv")))
}


files.ls <- list.files(path = count.input.dir, pattern = "csv", full.names = FALSE)
project_matrix_list <- NULL

for (f in files.ls) {
  temp_name = substring(f, 1, (nchar(f)-4))
  print(temp_name)
  temp_file = read.csv(file.path(count.input.dir, f), header = FALSE, sep = "\t")[, c(1,4)]
  colnames(temp_file) = c("ID", paste0(temp_name))
  temp_file = column_to_rownames(temp_file, "ID")
  project_matrix_list[[temp_name]] = temp_file
}

df <- project_matrix_list[[1]]
for (i in 2:length(project_matrix_list)) {
  print(names(project_matrix_list)[i])
  df <- merge(df, project_matrix_list[[i]], by = "row.names", all = TRUE)
  rownames(df) <- df$Row.names
  df <- df[ , !(names(df) %in% "Row.names")]  
}

reference = read.csv(file.path(count.input.dir, f), header = FALSE, sep = "\t")[,c(1,2)]
colnames(reference) <- c("ID", "Symbol")

df <- rownames_to_column(df, "ID")
df <- full_join(df, reference, by = "ID")
df <- df[-c(1:5),c(1,38,2:37)]

saveRDS(df, file.path(input.dir, "counts.RDS"))