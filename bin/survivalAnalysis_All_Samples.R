rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(futile.logger)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(stringr)
})


flog.threshold(DEBUG)


flog.debug("Set project directory")

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()
data.dir = file.path(proj.dir, "data", "TCGA", "survivalAnalysis_All")
dir.create(data.dir)


flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
FileNameData <- paste0(DataDirectory, platform,".rda")
file.type <- "results"


flog.debug("Get clinical data")

dataClin_COAD <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
dataClin_READ <- GDCquery_clinic(project = "TCGA-READ", type = "clinical") 

# join by common columns
common_col_names <- intersect(colnames(dataClin_READ), colnames(dataClin_COAD))
dataClin <- merge(dataClin_COAD, dataClin_READ, by = common_col_names, all = TRUE) 
which(duplicated(dataClin))
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)
levels <- dataClin$ajcc_pathologic_stage
dataClin$days_to_death <- as.numeric(dataClin$days_to_death)


flog.debug("Query GDC")

query <- GDCquery(project = CancerProject, 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = platform,
                  file.type = file.type,
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)

SampleNT <- TCGAquery_SampleTypes(getResults(query, cols = "cases"), "NT") 
SampleTP <- TCGAquery_SampleTypes(getResults(query, cols = "cases"), "TP") 


flog.debug("Prepare data and perform DEA")

UP = 0.6; DOWN = -0.6; FDR_cutoff = 0.05; reads = 50;
PreProc_cor.cut = 0.6; 
Norm_method = "gcContent"; 
Filt_method = "quantile"; 
Filt_qnt.cut = 0.25;
DEA_batch.factor = "Plate"; 
DEA_method = "glmLRT"


flog.debug("Query GDC")
queryDown <- GDCquery(
  project = CancerProject,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = platform,
  file.type = file.type,
  barcode = c(SampleTP, SampleNT),
  experimental.strategy = "RNA-Seq",
  legacy = TRUE)


flog.debug("Download samples")
tryCatch(GDCdownload(query = queryDown,
                     method = "api", 
                     files.per.chunk = 20,
                     directory = file.path(data.dir, "GDCdata")),
         error = function(e) GDCdownload(query = queryDown,
                                         method = "client", 
                                         files.per.chunk = 20,
                                         directory = file.path(data.dir, "GDCdata")))


flog.debug("Prepare GDC data")
FileNameData <- paste0("GDCdata", str_replace_all(platform, fixed(" "), ""),".rda")
dataPrep <- GDCprepare(
  query = queryDown, 
  save = TRUE,
  directory = file.path(data.dir, "GDCdata"), 
  save.filename = FileNameData)


flog.debug("Perform intensity correlation")
dataPreProc <- TCGAanalyze_Preprocessing(
  object = dataPrep, 
  cor.cut = PreProc_cor.cut)


flog.debug("Perform normalization")
dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPreProc,
  geneInfo = geneInfo,
  method = Norm_method)


flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm, 
  method = Filt_method, 
  qnt.cut =  Filt_qnt.cut)


flog.debug("Perform differential gene expression analysis")
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,SampleNT], 
  mat2 = dataFilt[,SampleTP],
  Cond1type = "Normal", 
  Cond2type = "Tumor",
  batch.factors = DEA_batch.factor,
  fdr.cut = 0.05, 
  logFC.cut = 0.6,
  method = DEA_method)


flog.debug("Perform survival analysis")

dataSurvival <- TCGAanalyze_SurvivalKM(
  clinical_patient = dataClin,
  dataGE = dataFilt,
  Genelist = c("VPS37A", "VPS37B", "VPS37C"),
  Survresult = FALSE,
  ThreshTop = 0.67,
  ThreshDown = 0.33,
  p.cut = 0.05, 
  group1 = SampleNT, 
  group2 = SampleTP) 


dataSurvival_all <- TCGAanalyze_SurvivalKM(
  clinical_patient = dataClin,
  dataGE = dataFilt,
  Genelist = rownames(dataDEGs),
  Survresult = FALSE,
  ThreshTop = 0.67,
  ThreshDown = 0.33,
  p.cut = 0.05, 
  group1 = SampleNT, 
  group2 = SampleTP) 

dataSurvival_all <- rownames_to_column(dataSurvival_all, "GeneSymbol")

saveRDS(dataSurvival_all, file.path(data.dir, "dataSurvival_all_Samples.RDS"))

sessionInfo()