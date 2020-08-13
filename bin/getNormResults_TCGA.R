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
  library(WriteXLS)
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


flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
FileNameData <- paste0(DataDirectory, platform,".rda")
file.type <- "normalized_results"

genes <- read.table(
  file.path("~/Desktop/HT projects/CRC_transcriptomics_proj/data/endocytic_genes.txt")) %>% 
  t() %>% c()
genes <- substring(genes, 1, (nchar(genes)-1))

dataClin <- readRDS(file.path(data.dir, "ClinData.RDS"))
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)
dataClin$ajcc_pathologic_t <- as.factor(dataClin$ajcc_pathologic_t)

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

DEGanalysis <- function(stages, title, 
                        reads = 50,
                        PreProc_cor.cut = 0.6, 
                        Filt_method = "quantile", 
                        Filt_qnt.cut = 0.25) {
  
  title = str_replace_all(title, fixed(" "), "")
  
  # Select only samples with clinical information
  flog.debug("Select only samples with clinical information")
  for(s in stages) {
    Stage_NT = filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleNT_Stage = filter(Stage_NT, Stage_NT$bcr_patient_barcode %in% Sample.Short)
    SampleNT_Stage = pmatch(SampleNT_Stage$bcr_patient_barcode, SampleNT)
    SampleNT_Stage = SampleNT[SampleNT_Stage]
    
    SampleNT_Stage_final = c(SampleNT_Stage_final, SampleNT_Stage)
    
    Stage_TP = filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleTP_Stage = filter(Stage_TP, Stage_TP$bcr_patient_barcode %in% Sample.Short)
    SampleTP_Stage = pmatch(SampleTP_Stage$bcr_patient_barcode, SampleTP)
    SampleTP_Stage = SampleTP[SampleTP_Stage]
    
    SampleTP_Stage_final = c(SampleTP_Stage_final, SampleTP_Stage)
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
  FileNameData = paste0("GDCdata", str_replace_all(platform, fixed(" "), ""),".rda")
  dataPrep = GDCprepare(
    query = queryDown, save = TRUE,
    directory = file.path(data.dir, "GDCdata"), 
    save.filename = FileNameData)
  saveRDS(dataPrep, file.path(data.dir, paste0("GDCprepare_", title, ".RDS")))
  # Samples preprocessing
  flog.debug("Perform intensity correlation")
  dataPreProc = TCGAanalyze_Preprocessing(
    object = dataPrep, 
    cor.cut = PreProc_cor.cut)
  saveRDS(dataPreProc, file.path(data.dir, paste0("dataPreProc_", title, ".RDS")))
  # Data filtering
  flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
  dataFilt = assign(paste0("dataFilt_", title), 
                     TCGAanalyze_Filtering(
                       tabDF = dataPreProc,
                       method = Filt_method, 
                       qnt.cut =  Filt_qnt.cut))
  dataFilt_transposed = t(dataFilt)
  dataFilt_transposed = as.data.frame(dataFilt_transposed)
  dataFilt_transposed = rownames_to_column(dataFilt_transposed, "barcode_id")
  # Selecting filtered data based on sample type
  F_SampleNT = TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "NT")
  F_SampleTP = TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
  ###  normal tissue
  dataFilt_transposed_NT = filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% F_SampleNT)
  dataFilt_transposed_NT = column_to_rownames(dataFilt_transposed_NT, "barcode_id")
  dataFilt_transposed_NT = t(dataFilt_transposed_NT)
  dataFilt_transposed_NT = as.data.frame(dataFilt_transposed_NT)
  temp_colnames = colnames(dataFilt_transposed_NT)
  temp_colnames = paste0("Colon_", temp_colnames)
  colnames(dataFilt_transposed_NT) = temp_colnames
  dataFilt_transposed_NT = rownames_to_column(dataFilt_transposed_NT, "Symbol")
  ### select only genes, which have at least x counts in all patients
  filter_NT = (rowSums(
    dataFilt_transposed_NT[,2:ncol(dataFilt_transposed_NT)]) >= ((ncol(dataFilt_transposed_NT)-1)*reads))
  dataFilt_transposed_NT$filter = filter_NT
  dataFilt_transposed_NT = filter(dataFilt_transposed_NT, dataFilt_transposed_NT$filter == TRUE)
  dataFilt_transposed_NT = dataFilt_transposed_NT[, -ncol(dataFilt_transposed_NT)]
  
  ### tumour tissue
  dataFilt_transposed_TP = filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% F_SampleTP)
  dataFilt_transposed_TP = column_to_rownames(dataFilt_transposed_TP, "barcode_id")
  dataFilt_transposed_TP = t(dataFilt_transposed_TP)
  dataFilt_transposed_TP = as.data.frame(dataFilt_transposed_TP)
  temp_colnames = colnames(dataFilt_transposed_TP)
  temp_colnames = paste0("Cancer_", temp_colnames)
  colnames(dataFilt_transposed_TP) = temp_colnames
  dataFilt_transposed_TP = rownames_to_column(dataFilt_transposed_TP, "Symbol")
  ### select only genes, which have at least x counts in all patients
  filter_TP = rowSums(
    dataFilt_transposed_TP[,2:ncol(dataFilt_transposed_TP)]) >= ((ncol(dataFilt_transposed_TP)-1)*reads)
  dataFilt_transposed_TP$filter = filter_TP
  dataFilt_transposed_TP = filter(dataFilt_transposed_TP, dataFilt_transposed_TP$filter == TRUE)
  dataFilt_transposed_TP = dataFilt_transposed_TP[, -ncol(dataFilt_transposed_TP)]
  
  # combine dataframes
  flog.debug("Prepare dataframe with filtered data for patients included in the analysis")
  patients = assign(paste0("dataFilt_patients_", title), 
                    merge(dataFilt_transposed_NT, 
                          dataFilt_transposed_TP, 
                          by = "Symbol"))
  patients = as.matrix(column_to_rownames(patients, "Symbol"))
  saveRDS(patients,
          file.path(data.dir, paste0("patients_", title,".RDS")))
  
  flog.debug("Select genes related to endocytic transport")
  n.counts = rownames_to_column(as.data.frame(patients), "Symbol")
  n.counts_subset =  filter(
    n.counts, 
    n.counts$Symbol %in% genes)
  assign(paste0("n.counts_", title), n.counts_subset)
  saveRDS(n.counts_subset, file.path(data.dir, paste0("n.counts_subset_", title, ".RDS")))
  
  genes_subset = n.counts_subset$Symbol
  viz_input = setNames(vector("list", length(genes_subset)), genes_subset)
  
  for(i in 1:nrow(n.counts_subset)) {
    temp_gene = n.counts_subset[i,1]
    temp_data = n.counts_subset[i,]
    temp_data = as.data.frame(t(temp_data))
    temp_data = rownames_to_column(temp_data, "metadata")
    temp_column = c("newID", temp_gene)  # as.character(unname(unlist(as.data.frame(t(temp_data[1,])))))
    colnames(temp_data) = temp_column
    temp_data = temp_data[-1,]
    rownames(temp_data) = NULL
    temp_data$Status <- factor(ifelse(
      test = str_detect(temp_data$newID, "Cancer"),
      yes = "Cancer",
      no = "Colon"),
      levels = c("Colon", "Cancer"))
    temp_data[,2] = as.numeric(as.character(temp_data[,2]))
    viz_input[[temp_gene]] = temp_data
  }
  assign(paste0("Viz_input_", title), viz_input)
  saveRDS(viz_input, file.path(data.dir, paste0("Viz_input_", title, ".RDS")))
  
  flog.debug("Perform statistical analysis of healthy vs. cancer samples")
  stats = setNames(vector("list", length(genes_subset)), genes_subset)
  for (i in 1:length(viz_input)) {
    temp_data = viz_input[[i]]
    test = pairwise.t.test(
      temp_data[,2], 
      temp_data[,3], 
      paired = TRUE,
      p.adj = "hochberg")
    stats[[i]] = test
  }
  
  stats_df <- do.call(rbind.data.frame, stats)[, 3, drop = FALSE]
  stats_df <- rownames_to_column(stats_df, "Symbol")
  stats_df_pval_subset <- filter(stats_df, stats_df$p.value < 0.05)
  assign(paste0("Stat_Analysis_", title), stats_df_pval_subset)
  saveRDS(stats_df_pval_subset, file.path(data.dir, paste0("Stat_Analysis_", title, ".RDS")))
}


flog.debug("Perform analysis for the early stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage I", "Stage IA", 
                       "Stage II", "Stage IIA", "Stage IIB", "Stage IIC"), 
            title = "Early stages")


flog.debug("Perform analysis for the late stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", 
                       "Stage IV", "Stage IVA", "Stage IVB"), 
            title = "Late stages")
