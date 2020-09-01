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
  library(biomaRt)
  library(WriteXLS)
})

### FUNCTIONS
`%!in%` <- negate(`%in%`)

### GLOBAL RESOURCES
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))


flog.threshold(DEBUG)

flog.debug("Set project directory")

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()
data.dir = file.path(proj.dir, "data", "TCGA", "RNASeq-Count-subset")
dir.create(data.dir, recursive = TRUE)


flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
WorkflowType = "HTSeq - Counts"


genes <- read.table(
  file.path("~/Desktop/HT projects/CRC_transcriptomics_proj/data/endocytic_genes.txt")) %>% 
  t() %>% c()
genes <- substring(genes, 1, (nchar(genes)-1))


dataClin_COAD <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical") 
dataClin_READ <- GDCquery_clinic(project = "TCGA-READ", type = "clinical") 
common_col_names <- intersect(colnames(dataClin_READ), colnames(dataClin_COAD))
dataClin <- merge(dataClin_COAD, dataClin_READ, by = common_col_names, all = TRUE) 
which(duplicated(dataClin))
dataClin <- dataClin[complete.cases(dataClin[, 15]), ]
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)

stats <- dataClin %>% dplyr::group_by(ajcc_pathologic_stage) %>% dplyr::summarise(n = n())
stats <- filter(stats, stats$ajcc_pathologic_stage %!in% NA)
stats <- stats$ajcc_pathologic_stage
stats <- as.character(stats)

dataClin <- filter(dataClin, dataClin$ajcc_pathologic_stage %in% stats)

saveRDS(dataClin, file.path(data.dir, paste0("dataClin.RDS")))


flog.debug("Query GDC")

query <- GDCquery(
  project = CancerProject, 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - Counts",
  legacy = FALSE) 

samplesDown <- query$results[[1]]$cases

MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(samplesDown, c("NT","TP"))
Sample.Short <- unique(substr(x = MatchedCoupledSampleTypes, start = 1, stop = 12))

SampleTP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")
SampleTP_pmatched <- pmatch(Sample.Short, SampleTP)
SampleTP_pmatched <- SampleTP_pmatched[!is.na(SampleTP_pmatched)]
SampleTP_final <- SampleTP[SampleTP_pmatched]

Sample.Short <- substr(x = SampleTP_final, start = 1, stop = 12)

SampleNT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
SampleNT_pmatched <- pmatch(Sample.Short, SampleNT)
SampleNT_pmatched <- SampleNT_pmatched[!is.na(SampleNT_pmatched)]
SampleNT_final <- SampleNT[SampleNT_pmatched]

### Load samples that were analyzed with protocol legacy = TRUE
Sample_Earlystages_Analyzed <- readRDS("~/Desktop/HT projects/TCGA_transcriptomics_proj/data_results/SampleNTEarlystages.RDS")
Sample_Latestages_Analyzed <- readRDS("~/Desktop/HT projects/TCGA_transcriptomics_proj/data_results/SampleNTLatestages.RDS")
Sample_Analyzed <- c(Sample_Earlystages_Analyzed, Sample_Latestages_Analyzed)
Sample_Analyzed <- substr(x = Sample_Analyzed, start = 1, stop = 12)
Sample_Analyzed <- paste(Sample_Analyzed, collapse = "|")


SampleTP_final <- data.frame(barcode = SampleTP_final)
SampleTP_final <- SampleTP_final %>%
  filter(!str_detect(barcode, Sample_Analyzed))


SampleNT_final <- data.frame(barcode = SampleNT_final)
SampleNT_final <- SampleNT_final %>%
  filter(!str_detect(barcode, Sample_Analyzed))

Sample.Short <- substr(x = SampleNT_final$barcode, start = 1, stop = 12)

dataClin <- filter(
  dataClin,
  dataClin$bcr_patient_barcode %in% Sample.Short)

dataClin %>% dplyr::group_by(ajcc_pathologic_stage) %>% dplyr::summarise(n = n())

Sample.Short_4_filter <- dataClin$submitter_id
Sample.Short_4_filter <- paste0(Sample.Short, collapse = "|")



flog.debug("Download, normalize, filter and perform differential gene expression analysis")

DEGanalysis <- function(stages, title,
                        UP = 0.6, DOWN = -0.6, FDR_cutoff = 0.05,
                        reads = 5,
                        PreProc_cor.cut = 0.6, 
                        Filt_method = "quantile", 
                        Filt_qnt.cut = 0.25,
                        DEA_batch.factor = "Plate",
                        DEA_method = "glmLRT") {
  
  title = str_replace_all(title, fixed(" "), "")
  
  # Select only samples with clinical information
  flog.debug("Select only samples with clinical information")
  for(s in stages) {
    Stage_NT = filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleNT_Stage = filter(Stage_NT, Stage_NT$bcr_patient_barcode %in% Sample.Short)
    SampleNT_Stage = pmatch(SampleNT_Stage$bcr_patient_barcode, SampleNT)
    SampleNT_Stage = SampleNT[SampleNT_Stage]
    
    SampleNT_Stage_final = c(SampleNT_Stage_final, SampleNT_Stage)
    
    saveRDS(SampleNT_Stage_final, file.path(data.dir, paste0("SampleNT", title, ".RDS")))
    
    Stage_TP = filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleTP_Stage = filter(Stage_TP, Stage_TP$bcr_patient_barcode %in% Sample.Short)
    SampleTP_Stage = pmatch(SampleTP_Stage$bcr_patient_barcode, SampleTP)
    SampleTP_Stage = SampleTP[SampleTP_Stage]
    
    SampleTP_Stage_final = c(SampleTP_Stage_final, SampleTP_Stage)
    
    saveRDS(SampleTP_Stage_final, file.path(data.dir, paste0("SampleTP", title, ".RDS")))
  }
  
  # Query platform Illumina HiSeq to download samples
  flog.debug("Query GDC")
  assign("queryDown", 
         GDCquery(
           project = CancerProject, 
           data.category = "Transcriptome Profiling",
           data.type = "Gene Expression Quantification", 
           workflow.type = "HTSeq - Counts",
           legacy = FALSE))
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
  FileNameData = paste0("GDCdata", str_replace_all(WorkflowType, fixed(" "), ""),".rda")
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
  ###  normal tissue
  dataFilt_transposed_NT = filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% SampleNT_Stage_final)
  dataFilt_transposed_NT = column_to_rownames(dataFilt_transposed_NT, "barcode_id")
  dataFilt_transposed_NT = t(dataFilt_transposed_NT)
  dataFilt_transposed_NT = as.data.frame(dataFilt_transposed_NT)
  temp_colnames = colnames(dataFilt_transposed_NT)
  temp_colnames = paste0("Colon_", temp_colnames)
  colnames(dataFilt_transposed_NT) = temp_colnames
  dataFilt_transposed_NT = rownames_to_column(dataFilt_transposed_NT, "ENSEMBLid")
  ### select only genes, which have at least x counts in all patients
  filter_NT = (rowSums(
    dataFilt_transposed_NT[,2:ncol(dataFilt_transposed_NT)]) >= ((ncol(dataFilt_transposed_NT)-1)*reads))
  dataFilt_transposed_NT$filter = filter_NT
  dataFilt_transposed_NT = filter(dataFilt_transposed_NT, dataFilt_transposed_NT$filter == TRUE)
  dataFilt_transposed_NT = dataFilt_transposed_NT[, -ncol(dataFilt_transposed_NT)]
  
  ### tumour tissue
  dataFilt_transposed_TP = filter(
    dataFilt_transposed, 
    dataFilt_transposed$barcode_id %in% SampleTP_Stage_final)
  dataFilt_transposed_TP = column_to_rownames(dataFilt_transposed_TP, "barcode_id")
  dataFilt_transposed_TP = t(dataFilt_transposed_TP)
  dataFilt_transposed_TP = as.data.frame(dataFilt_transposed_TP)
  temp_colnames = colnames(dataFilt_transposed_TP)
  temp_colnames = paste0("Cancer_", temp_colnames)
  colnames(dataFilt_transposed_TP) = temp_colnames
  dataFilt_transposed_TP = rownames_to_column(dataFilt_transposed_TP, "ENSEMBLid")
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
                          by = "ENSEMBLid"))
  patients = as.matrix(column_to_rownames(patients, "ENSEMBLid"))
  saveRDS(patients,
          file.path(data.dir, paste0("patients_", title,".RDS")))
  
  flog.debug("Perform differential gene expression analysis")
  dataDEGs <- TCGAanalyze_DEA(mat1 = patients[,paste0("Colon_", SampleNT_Stage_final)], 
                              mat2 = patients[,paste0("Cancer_", SampleTP_Stage_final)], 
                              Cond1type = "Normal", 
                              Cond2type = "Tumor",
                              batch.factors = DEA_batch.factor,
                              fdr.cut = 1, 
                              logFC.cut = 0,
                              method = DEA_method)
  
  dataDEGs <- tibble::rownames_to_column(dataDEGs, var = "ENSEMBLid")
  dataDEGs <- mutate_at(dataDEGs, vars(-ENSEMBLid), funs(as.numeric(.)))
  dataDEGs$Threshold <- with(dataDEGs, ifelse(dataDEGs$logFC >= UP & 
                                                dataDEGs$FDR < FDR_cutoff, 
                                              "Upregulated",
                                              ifelse(dataDEGs$logFC <= DOWN & 
                                                       dataDEGs$FDR < FDR_cutoff, 
                                                     "Downregulated", "Not significant")))
  dataDEGs <<- assign(paste0("dataDEGs_", title), dataDEGs)
  saveRDS(dataDEGs,
          file.path(data.dir, paste0("dataDEGs_", title,".RDS")))
  
  flog.debug("Convert ENSEMBLid to Gene Symbols")
  ENSEMBLid = dataDEGs$ENSEMBLid
  
  Symbols = getBM(
    filters = "ensembl_gene_id", 
    attributes = c("ensembl_gene_id","hgnc_symbol"), 
    values = ENSEMBLid, 
    mart = mart)
  
  Symbols = dplyr::rename(
    Symbols, 
    ENSEMBLid = ensembl_gene_id,
    Symbol = hgnc_symbol)
  
  dataDEGs = inner_join(dataDEGs, Symbols, by = "ENSEMBLid") 
  
  saveRDS(dataDEGs, file.path(data.dir, paste0("dataDEGs_", title, ".RDS")))
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
