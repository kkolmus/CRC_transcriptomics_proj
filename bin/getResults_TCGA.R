rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
  library(WriteXLS)
})


flog.threshold(DEBUG)

flog.debug("Set working environment")

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()
data.dir <- file.path(proj.dir, "data", "TCGA", "Illumina_HiSeq_Results")


flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
FileNameData <- paste0(DataDirectory, platform,".rda")
file.type <- "results"

dataClin <- readRDS(file.path(data.dir, "ClinData_matchedSamples.RDS"))
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)


flog.debug("Query GDC")

query <- GDCquery(
  project = CancerProject, 
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = platform,
  file.type = file.type,
  experimental.strategy = "RNA-Seq",
  legacy = TRUE)

samplesDown <- query$results[[1]]$cases

MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(samplesDown, c("NT","TP"))
SampleNT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
SampleTP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")
Sample.Short <- unique(substr(x = MatchedCoupledSampleTypes, start = 1, stop = 12))
rm(samplesDown)


flog.debug("Download, normalize, filter and perform differential gene expression analysis")

DEGanalysis <- function(
  stages, 
  title,
  UP = 0.6, 
  DOWN = -0.6, 
  confidence.cutoff = 0.05, 
  reads = 50,
  PreProc_cor.cut = 0.6, 
  Norm_method = "gcContent", 
  Filt_method = "quantile", 
  Filt_qnt.cut = 0.25) {
  
  title = str_replace_all(title, fixed(" "), "")
  
  # Select only samples with clinical information
  flog.debug("Select only samples with clinical information")
  for(s in stages) {
    Stage_NT <- filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleNT_Stage <- filter(Stage_NT, Stage_NT$bcr_patient_barcode %in% Sample.Short)
    SampleNT_Stage <- pmatch(SampleNT_Stage$bcr_patient_barcode, SampleNT)
    SampleNT_Stage <- SampleNT[SampleNT_Stage]
    
    SampleNT_Stage_final <- c(SampleNT_Stage_final, SampleNT_Stage)
    
    saveRDS(SampleNT_Stage_final, file.path(data.dir, paste0("SampleNT", title, ".RDS")))
    
    Stage_TP <- filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleTP_Stage <- filter(Stage_TP, Stage_TP$bcr_patient_barcode %in% Sample.Short)
    SampleTP_Stage <- pmatch(SampleTP_Stage$bcr_patient_barcode, SampleTP)
    SampleTP_Stage <- SampleTP[SampleTP_Stage]
    
    SampleTP_Stage_final <- c(SampleTP_Stage_final, SampleTP_Stage)
    
    saveRDS(SampleTP_Stage_final, file.path(data.dir, paste0("SampleTP", title, ".RDS")))
  }
  
  # Query platform Illumina HiSeq to download samples
  flog.debug("Query GDC")
  assign("queryDown", 
         GDCquery(
           project = CancerProject,
           data.category = "Gene expression",
           data.type = "Gene expression quantification",
           platform = platform,
           file.type = file.type,
           barcode = c(SampleTP_Stage_final, SampleNT_Stage_final),
           experimental.strategy = "RNA-Seq",
           legacy = TRUE))
  saveRDS(queryDown, file.path(data.dir, paste0("GDCquery_", title, ".RDS")))
  # Download samples
  flog.debug("Download samples")
  tryCatch(GDCdownload(
    query = queryDown,
    method = "api", 
    files.per.chunk = 20,
    directory = file.path(data.dir, "GDCdata")),
    error = function(e) GDCdownload(
      query = queryDown,
      method = "client", 
      files.per.chunk = 20,
      directory = file.path(data.dir, "GDCdata")))
  # Prepare samples for analysis
  flog.debug("Prepare GDC data")
  FileNameData <- paste0("GDCdata", str_replace_all(platform, fixed(" "), ""),".rda")
  dataPrep <- GDCprepare(
    query = queryDown, 
    save = TRUE,
    directory = file.path(data.dir, "GDCdata"), 
    save.filename = FileNameData)
  saveRDS(dataPrep, file.path(data.dir, paste0("GDCprepare_", title, ".RDS")))
  # Samples preprocessing
  flog.debug("Perform intensity correlation")
  dataPreProc <- TCGAanalyze_Preprocessing(
    object = dataPrep, 
    cor.cut = PreProc_cor.cut)
  saveRDS(dataPreProc, file.path(data.dir, paste0("dataPreProc_", title, ".RDS")))
  # Samples normalisation
  flog.debug("Perform normalization")
  dataNorm <- TCGAanalyze_Normalization(
    tabDF = dataPreProc,
    geneInfo = geneInfo,
    method = Norm_method)
  saveRDS(dataNorm, file.path(data.dir, paste0("dataNorm_", title, ".RDS")))
  # Data filtering
  flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
  dataFilt <- assign(
    paste0("dataFilt_", title), 
    TCGAanalyze_Filtering(
      tabDF = dataNorm, 
      method = Filt_method, 
      qnt.cut =  Filt_qnt.cut))
  dataFilt_transposed <- t(dataFilt)
  dataFilt_transposed <- as.data.frame(dataFilt_transposed)
  dataFilt_transposed <- rownames_to_column(dataFilt_transposed, "barcode_id")
  # Selecting filtered data based on sample type
  F_SampleNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "NT")
  F_SampleTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
  # normal tissue
  dataFilt_transposed_NT <- filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% F_SampleNT)
  dataFilt_transposed_NT <- column_to_rownames(dataFilt_transposed_NT, "barcode_id")
  dataFilt_transposed_NT <- t(dataFilt_transposed_NT)
  dataFilt_transposed_NT <- as.data.frame(dataFilt_transposed_NT)
  dataFilt_transposed_NT <- rownames_to_column(dataFilt_transposed_NT, "Symbol")
  # select only genes, which have at least 10 counts in all patients
  filter_NT <- (rowSums(dataFilt_transposed_NT[,2:ncol(dataFilt_transposed_NT)]) >= ((ncol(dataFilt_transposed_NT)-1)*reads))
  dataFilt_transposed_NT$filter <- filter_NT
  dataFilt_transposed_NT <- filter(dataFilt_transposed_NT, dataFilt_transposed_NT$filter == TRUE)
  dataFilt_transposed_NT <- dataFilt_transposed_NT[, -ncol(dataFilt_transposed_NT)]
  
  # tumour tissue
  dataFilt_transposed_TP <- filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% F_SampleTP)
  dataFilt_transposed_TP <- column_to_rownames(dataFilt_transposed_TP, "barcode_id")
  dataFilt_transposed_TP <- t(dataFilt_transposed_TP)
  dataFilt_transposed_TP <- as.data.frame(dataFilt_transposed_TP)
  dataFilt_transposed_TP <- rownames_to_column(dataFilt_transposed_TP, "Symbol")
  # select only genes, which have at least 10 counts in all patients
  filter_TP <- rowSums(dataFilt_transposed_TP[,2:ncol(dataFilt_transposed_TP)]) >= ((ncol(dataFilt_transposed_TP)-1)*reads)
  dataFilt_transposed_TP$filter <- filter_TP
  dataFilt_transposed_TP <- filter(dataFilt_transposed_TP, dataFilt_transposed_TP$filter == TRUE)
  dataFilt_transposed_TP <- dataFilt_transposed_TP[, -ncol(dataFilt_transposed_TP)]
  
  # combine dataframes
  flog.debug("Prepare dataframe with filtered data for patients included in the analysis")
  patients <- assign(
    paste0("dataFilt_patients_", title), 
    merge(
      dataFilt_transposed_NT, 
      dataFilt_transposed_TP, by = "Symbol"))
  patients <- as.matrix(column_to_rownames(patients, "Symbol"))
  saveRDS(patients, file.path(data.dir, paste0("patients_", title,".RDS")))
  # Differential gene expression analysis 
  flog.debug("Perform differential gene expression analysis")
  # prepare count matrix
  mat1 = patients[,F_SampleTP]
  new_colnames = paste0(colnames(mat1), "_", "cancer")
  colnames(mat1) = new_colnames
  dim(mat1)
  mat2 = patients[,F_SampleNT]
  new_colnames = paste0(colnames(mat2), "_", "colon")
  colnames(mat2) = new_colnames
  dim(mat2)
  temp_cts = cbind(mat1, mat2)
  dim(temp_cts)
  rm(new_colnames, mat1, mat2)
  # prepare pheno data
  temp_pdata = data.frame(
    sampleID = colnames(temp_cts))
  temp_pdata$tissue_type = factor(ifelse(
    test = grepl(pattern = "cancer", temp_pdata$sampleID),
    yes = "cancer",
    no = "colon"
  ), levels = c("colon", "cancer"))
  temp_pdata$plate = factor(substring(temp_pdata$sampleID, 22, 25))
  rownames(temp_pdata) = temp_pdata$sampleID
  
  dds = DESeqDataSetFromMatrix(
    countData = temp_cts,
    colData = temp_pdata,
    design = ~ tissue_type + plate)

  dds = DESeq(dds)
  res = results(dds)
  res = results(dds, name = "tissue_type_cancer_vs_colon")
  res_gene = res@rownames
  res = as.data.frame(res@listData)
  rownames(res) = res_gene
  
  dataDEAs = tibble::rownames_to_column(res, var = "Symbol")
  dataDEAs = mutate_at(dataDEAs, vars(-Symbol), funs(as.numeric(.)))
  dataDEAs$Threshold <- with(
    dataDEAs, 
    ifelse(
      test = dataDEAs$log2FoldChange >= UP & dataDEAs$padj < confidence.cutoff, 
      yes = "Upregulated",
      no = ifelse(
        test = dataDEAs$log2FoldChange <= DOWN & dataDEAs$padj < confidence.cutoff, 
        yes = "Downregulated", 
        no = "Not significant")))
  
  temp_cts = as.data.frame(temp_cts)
  temp_cts = rownames_to_column(temp_cts, "Symbol")
  
  dataDEAs = left_join(dataDEAs, temp_cts, by = "Symbol")
  
  dataDEAs <<- assign(paste0("dataDEAs_", title), dataDEAs)
  
  saveRDS(dataDEAs, file.path(data.dir, paste0("dataDEA_", title,".RDS")))
}


flog.debug("Perform analysis for the early stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(
  stages = c("Stage I", "Stage IA",  "Stage II", "Stage IIA", "Stage IIB", "Stage IIC"), 
  title = "Early stages")


flog.debug("Perform analysis for the late stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(
  stages = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA", "Stage IVB"), 
  title = "Late stages")


flog.debug("Combine datasets")

dataDEA_Earlystages <- readRDS(file.path(data.dir, "dataDEA_Earlystages.RDS"))
dataDEA_Latestages <- readRDS(file.path(data.dir, "dataDEA_Latestages.RDS"))

dataDEA_Earlystages <- dataDEA_Earlystages[, c(1,3,7,8)]
dataDEA_Earlystages <- filter(dataDEA_Earlystages, dataDEA_Earlystages$Threshold != "Not significant")

dataDEA_Latestages <- dataDEA_Latestages[, c(1,3,7,8)]
dataDEA_Latestages <- filter(dataDEA_Latestages, dataDEA_Latestages$Threshold != "Not significant")

dataDEA <- full_join(dataDEA_Earlystages, dataDEA_Latestages, by = "Symbol",
                     suffix = c("_earlyStage", "_advancedStages"))
dataDEA[is.na(dataDEA)] <- ""

setwd("~/Desktop/HT projects/CRC_transcriptomics_proj/data/TCGA/Illumina_HiSeq_Results")
WriteXLS(dataDEA, "dataDEA.xlsx")