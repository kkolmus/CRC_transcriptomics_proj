rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(futile.logger)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(tidyverse)
})

flog.threshold(DEBUG)


flog.debug("Set project directory")

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()
data.dir = file.path(proj.dir, "data", "TCGA", "survivalAnalysis")
dir.create(data.dir, recursive = TRUE)


flog.debug("Get clinical data")

dataClin_COAD <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
# 459 cases
dataClin_READ <- GDCquery_clinic(project = "TCGA-READ", type = "clinical") 
# 171 cases

# join by common columns
common_col_names <- intersect(colnames(dataClin_READ), colnames(dataClin_COAD))
dataClin <- merge(dataClin_COAD, dataClin_READ, by = common_col_names, all = TRUE) 
which(duplicated(dataClin))
# 630 cases

# remove columns with only missing data and select data for downstream analysis
dataClin <- dataClin[, colSums(is.na(dataClin)) != nrow(dataClin)] %>% 
  filter(dataClin$tumor_stage != "not reported")
# 609 cases

dataClin %>% group_by(dataClin$tumor_stage) %>% summarise(n = n())
# `summarise()` ungrouping output (override with `.groups` argument)
# # A tibble: 13 x 2
#   `dataClin$tumor_stage`      n
#   <chr>                   <int>
# 1 stage i                   108
# 2 stage ia                    1
# 3 stage ii                   38
# 4 stage iia                 177
# 5 stage iib                  12
# 6 stage iic                   2
# 7 stage iii                  25
# 8 stage iiia                 15
# 9 stage iiib                 86
# 10 stage iiic                55
# 11 stage iv                  64
# 12 stage iva                 24
# 13 stage ivb                  2

dataClin %>% group_by(year_of_diagnosis) %>% summarise(n = n())
# `summarise()` ungrouping output (override with `.groups` argument)
# # A tibble: 16 x 2
#    year_of_diagnosis     n
#                <int> <int>
# 1               1998     1
# 2               1999     4
# 3               2000    13
# 4               2001    14
# 5               2002    13
# 6               2003     8
# 7               2004    18
# 8               2005    26
# 9               2006    34
# 10              2007    73
# 11              2008    68
# 12              2009   117
# 13              2010   120
# 14              2011    77
# 15              2012    11
# 16              2013    12


dataClin %>% group_by(vital_status) %>% summarise(n = n())
# `summarise()` ungrouping output (override with `.groups` argument)
# # A tibble: 2 x 2
#   vital_status     n
#   <chr>        <int>
# 1 Alive          488
# 2 Dead           121



dataClin %>% group_by(year_of_death) %>% summarise(n = n())
# `summarise()` ungrouping output (override with `.groups` argument)
# # A tibble: 15 x 2
#     year_of_death     n
#             <int> <int>
# 1           1998     1
# 2           1999     1
# 3           2001     2
# 4           2002     7
# 5           2003     5
# 6           2004     5
# 7           2005     2
# 8           2006     5
# 9           2007     7
# 10          2008     6
# 11          2009    15
# 12          2010     2
# 13          2011     1
# 14          2012     2
# 15            NA   548


