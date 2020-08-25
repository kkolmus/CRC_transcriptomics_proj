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
data.dir = file.path(proj.dir, "data", "TCGA", "RNASeq")
dir.create(data.dir, recursive = TRUE)


flog.debug("List available TCGA projects")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects


flog.debug("Load information required for the analysis")

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
file.type <- "normalized_results"
FileNameData <- paste0(DataDirectory, platform,".rda")


flog.debug("Query platform Illumina HiSeq for list of barcodes")

query <- GDCquery(
  project = CancerProject, 
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = platform,
  file.type = file.type,
  experimental.strategy = "RNA-Seq",
  legacy = TRUE)

samplesDown <- query$results[[1]]$cases


flog.debug("Define sample type")

SampleTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
saveRDS(SampleTP, file.path(data.dir, "SampleTP.RDS"))
SampleNT <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")
saveRDS(SampleNT, file.path(data.dir, "SampleNT.RDS"))


flog.debug("Extract short barcode - project(4)-TSS(2)-participant(4)")

SampleTP_short <- substr(x = SampleTP, start = 1, stop = 12)
SampleNT_short <- substr(x = SampleNT, start = 1, stop = 12)


flog.debug("Get clinical data")

dataClin_COAD <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical") 
dataClin_READ <- GDCquery_clinic(project = "TCGA-READ", type = "clinical") 
# join by common columns
common_col_names <- intersect(colnames(dataClin_READ), colnames(dataClin_COAD))
dataClin <- merge(dataClin_COAD, dataClin_READ, by = common_col_names, all = TRUE) 
which(duplicated(dataClin))

# remove columns with only missing data and select data for downstream analysis
dataClin <- dataClin[, colSums(is.na(dataClin)) != nrow(dataClin)] %>%
  select(disease, bcr_patient_barcode, primary_diagnosis, tissue_or_organ_of_origin, 
         site_of_resection_or_biopsy, prior_treatment, prior_malignancy, 
         ajcc_pathologic_stage, ajcc_pathologic_t,
         gender, vital_status, race, ethnicity, bmi, year_of_birth,
         year_of_diagnosis, year_of_death, age_at_diagnosis, days_to_birth, days_to_death, 
         days_to_last_follow_up, 
         treatments_pharmaceutical_treatment_or_therapy,
         treatments_radiation_treatment_or_therapy)

saveRDS(dataClin, file.path(data.dir, "ClinData.RDS"))

# tumor_stage = ajcc_pathologic_stage


flog.debug("Prepare summary of data on staging")

new_names <- c("Staging", "Number of patients")

dataClin_TP <- filter(dataClin, dataClin$bcr_patient_barcode %in% SampleTP_short) 
SampleTP_short <- dataClin_TP$bcr_patient_barcode
table_TP <- dataClin_TP %>% group_by(ajcc_pathologic_stage) %>% summarise(n())
colnames(table_TP) <- new_names

# normal tissue

dataClin_NT <- filter(dataClin, dataClin$bcr_patient_barcode %in% SampleNT_short) 
SampleNT_short <- dataClin_NT$bcr_patient_barcode
table_NT <- dataClin_NT %>% group_by(ajcc_pathologic_stage) %>% summarise(n())
colnames(table_NT) <- new_names

# identification of stagings with at least 3 patients

table_TP <- filter(table_TP, table_TP$Staging != "NA", table_TP$`Number of patients` >= 3)

table_NT <- filter(table_NT, table_NT$Staging != "NA", table_NT$`Number of patients` >= 3)

# combining tables

stagings <- full_join(table_TP, table_NT, by = "Staging", suffix = c("_TP", "_NT"))
stagings[is.na(stagings)] <- 0

saveRDS(stagings, file.path(data.dir, "Staging_summary.RDS"))


flog.debug("Prepare data subset based on patient's mortality")

MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(samplesDown, c("NT","TP"))
SamplesMatched_NT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
saveRDS(SamplesMatched_NT, file.path(data.dir, "matchedSamples_NT.RDS"))

SamplesMatched_TP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")
saveRDS(SamplesMatched_TP, file.path(data.dir, "matchedSamples_TP.RDS"))

SampleMatched_short <- substr(x = SamplesMatched_NT, start = 1, stop = 12)
saveRDS(SampleMatched_short, file.path(data.dir, "matchedSamples_short.RDS"))

dataClin_Matched <- dplyr::filter(dataClin, dataClin$bcr_patient_barcode %in% SampleMatched_short)
saveRDS(dataClin_Matched, file.path(data.dir, "ClinData_matchedSamples.RDS"))