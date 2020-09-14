rm(list = ls())

### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(scran)
  library(tidyverse)
  library(cowplot)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(VennDiagram)
  library(biomaRt)
})

### PARAMETERS ###

UP = 0.6
DOWN = -0.6
confidence.cutoff = 0.05

### FUNCTIONS ###

listMarts()
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

dataset = patients_early

TPM_calculator <- function(dataset = DF) {
  IDs <- dataset$Symbol
  length = getBM(attributes = c("hgnc_symbol", "transcript_length"), 
                 filters = "hgnc_symbol", values = IDs, mart = ensembl)
  length = length[! duplicated(length$hgnc_symbol),]
  length = dplyr::select(length, Symbol = hgnc_symbol, transcript_length)
  df_temp = merge(dataset, length, by = c("Symbol"))
  df_temp = df_temp[, c(1, ncol(df_temp), 2:(ncol(df_temp)-1))]
  # normalize for gene length, TPM1 = reads per kilobase = RPK
  TPM1 <- mutate_at(df_temp, vars(-c(Symbol)), funs(./df_temp$transcript_length))[, -2]
  # normalize for sequencing depth, TPM2 = scaling factors
  TPM2 <- unname(unlist(apply(TPM1[, c(2:ncol(TPM1))], 2, sum)/10^(6)))
  # transcripts per million
  TPM_temp <- sweep(TPM1[,-1], MARGIN = 2, TPM2, FUN = "/")
  TPM_names <- TPM1[,1, drop = FALSE]
  TPM <<- cbind(TPM_names, TPM_temp) # TPM
  # scaled samples
  Av <- apply(TPM[, c(2:ncol(TPM))], 1, mean)
  StDev <- apply(TPM[, c(2:ncol(TPM))], 1, sd)
  TPM_HM <<- mutate_at(TPM, vars(-Symbol), funs(((.)-Av)/StDev)) # scaled
}


### MAIN ###

flog.threshold(DEBUG)

flog.debug("Set working environment")

data.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj/data/statistics"
dir.create(data.dir)
viz.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj/data/vizualization"
dir.create(viz.dir)

TCGA_dataDEA_Early <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/TCGA/Illumina_HiSeq_Results/dataDEA_Earlystages.RDS")
TCGA_dataDEA_Early <- TCGA_dataDEA_Early[,c(1,3,7,8)]
TCGA_dataDEA_Early <- TCGA_dataDEA_Early[complete.cases(TCGA_dataDEA_Early),]

GES50760_dataDEA_Late <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/input/GSE50760/DEA.RDS")
GES50760_dataDEA_Late <- GES50760_dataDEA_Late[,c(2,3,4)]
GES50760_dataDEA_Late <- GES50760_dataDEA_Late[complete.cases(GES50760_dataDEA_Late),]

TCGA_dataDEA_Late <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/TCGA/Illumina_HiSeq_Results/dataDEA_Latestages.RDS")
TCGA_dataDEA_Late <- TCGA_dataDEA_Late[,c(1,3,7)]
TCGA_dataDEA_Late <- TCGA_dataDEA_Late[complete.cases(TCGA_dataDEA_Late),]


flog.debug("Combine project")

DEA_Early <- TCGA_dataDEA_Early
# 13292
DEA_Late <- inner_join(
  GES50760_dataDEA_Late, 
  TCGA_dataDEA_Late, 
  by = "Symbol", 
  suffix = c("GSE", "TCGA"))

rm(GES50760_dataDEA_Late, TCGA_dataDEA_Early, TCGA_dataDEA_Late)


flog.debug("Meta-analysis for advanced stages")

DEA_Late$log2FoldChange <- (DEA_Late$log2FoldChangeGSE + DEA_Late$log2FoldChangeTCGA)/2
DEA_Late$padj <- combinePValues(DEA_Late$padjGSE, DEA_Late$padjTCGA, method = "z")
DEA_Late$Threshold <- with(
  DEA_Late, 
  ifelse(
    test = DEA_Late$log2FoldChange >= UP & DEA_Late$padj < confidence.cutoff, 
    yes = "Upregulated",
    no = ifelse(
      test = DEA_Late$log2FoldChange <= DOWN & DEA_Late$padj < confidence.cutoff, 
      yes = "Downregulated", 
      no = "Not significant")))
DEA_Late <- DEA_Late[, -c(2:5)]
# 11086

flog.debug("Selection")

genes <- read.table(
  file.path("~/Desktop/HT projects/CRC_transcriptomics_proj/data/endocytic_genes.txt")) %>% 
  t() %>% c()
genes <- substring(genes, 1, (nchar(genes)-1))

DEA_Early_Filtered <- filter(DEA_Early, DEA_Early$Symbol %in% genes)
# 401 genes
DEA_Late_Filtered <- filter(DEA_Late, DEA_Late$Symbol %in% genes)
# 388 genes


flog.debug("Descriptive statistics for differentially expressed genes")

flog.debug("Differentially expressed genes in early stages of cancer progression")

early_up = filter(DEA_Early_Filtered, DEA_Early_Filtered$Threshold == "Upregulated")
saveRDS(early_up, file.path(data.dir, "stage_early_up.RDS"))
# 34

early_down = filter(DEA_Early_Filtered, DEA_Early_Filtered$Threshold == "Downregulated")
saveRDS(early_down, file.path(data.dir, "stage_early_down.RDS"))
# 66

early = rbind(early_up, early_down)
saveRDS(early, file.path(data.dir, "stage_early.RDS"))
# 100

late_up = filter(DEA_Late_Filtered, DEA_Late_Filtered$Threshold == "Upregulated")
saveRDS(late_up, file.path(data.dir, "stage_late_up.RDS"))
# 31

late_down = filter(DEA_Late_Filtered, DEA_Late_Filtered$Threshold == "Downregulated")
saveRDS(late_down, file.path(data.dir, "stage_late_down.RDS"))
# 38

late = rbind(late_up, late_down)
saveRDS(late, file.path(data.dir, "stage_late.RDS"))
# 69


flog.debug("Unique genes and genes at the intersection between stages")

intersect <- intersect(early$Symbol, late$Symbol)
cat("number of common genes: ", length(intersect))
# number of common genes :  57
saveRDS(intersect, file.path(data.dir, "intersect.RDS"))

unique_4_early <- setdiff(early$Symbol, intersect)
cat("number of unique genes in early stages: ", length(unique_4_early))
# number of unique genes in early stages:  43
saveRDS(unique_4_early, file.path(data.dir, "unique_4_early.RDS"))

unique_4_late <- setdiff(late$Symbol, intersect)
cat("number of unique genes in advanced stages: ", length(unique_4_late))
# number of unique genes in late stages:  12
saveRDS(unique_4_late, file.path(data.dir, "unique_4_late.RDS"))

rm(unique_4_late, unique_4_early, intersect, early, early_down, early_up, late, late_up, late_down)


flog.debug("Volcano plots")

GOI <- c("TSG101", "VPS28", "VPS37A", "VPS37B", "VPS37C", "VPS37D")

# early stages

dataDEGs <- DEA_Early
dataDEGs$Threshold <- gsub(pattern = "Not significant", "No change", dataDEGs$Threshold)
dataDEGs$Threshold <- gsub(pattern = "Downregulated", "Decreased", dataDEGs$Threshold)
dataDEGs$Threshold <- gsub(pattern = "Upregulated", "Increased", dataDEGs$Threshold)

GOI_filtered <- filter(dataDEGs, dataDEGs$Symbol %in% GOI)

vp <- ggplot(data = dataDEGs,
             aes(x = log2FoldChange,
                 y = -log10(padj),
                 colour = Threshold)) +
  scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
  geom_point(alpha = 0.4, size = 1.0) + xlim(c(-7.5, 7.5)) + ylim(c(0, 25)) +
  labs(color = "Expression: ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 14)) +
  geom_point(data = GOI_filtered, colour = "black") +
  # ggtitle(paste0("Early stages of carcinogenesis - Cancer vs. Colon")) +
  xlab("log2FoldChange") + ylab("-log10(p-value)") +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 7)))

vp + 
  geom_text_repel(data = GOI_filtered,
                  aes(label = GOI_filtered$Symbol), colour = "black", size = 5)



# late stages

dataDEGs <- DEA_Late
dataDEGs$Threshold <- gsub(pattern = "Not significant", "No change", dataDEGs$Threshold)
dataDEGs$Threshold <- gsub(pattern = "Downregulated", "Decreased", dataDEGs$Threshold)
dataDEGs$Threshold <- gsub(pattern = "Upregulated", "Increased", dataDEGs$Threshold)

GOI_filtered <- filter(dataDEGs, dataDEGs$Symbol %in% GOI)

vp <- ggplot(data = dataDEGs,
             aes(x = log2FoldChange,
                 y = -log10(padj),
                 colour = Threshold)) +
  scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
  geom_point(alpha = 0.4, size = 1.0) + xlim(c(-7.5, 7.5)) + ylim(c(0, 25)) +
  labs(color = "Expression: ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 14)) +
  geom_point(data = GOI_filtered, colour = "black") +
  # ggtitle(paste0("Late stages of carcinogenesis - Cancer vs. Colon")) +
  xlab("log2FoldChange") + ylab("-log10(p-value)") +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 7)))

vp + 
  geom_text_repel(data = GOI_filtered,
                  aes(label = GOI_filtered$Symbol), colour = "black", size = 5)

rm(vp, GOI_filtered, GOI, dataDEGs)


flog.debug("Heatmap")

patients_GSE50760 <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/input/GSE50760/DEA.RDS") %>%
  as.data.frame() %>%
  dplyr::select(-ID, -log2FoldChange, -padj) 

patients_lateStages <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/TCGA/Illumina_HiSeq_Results/patients_Latestages.RDS") %>%
  as.data.frame() %>%
  rownames_to_column("Symbol")

### MODIFY COLUMNS OF COUNT MATRIX FOR LATE STAGES

temp_colnames <- colnames(patients_lateStages)
new_temp_colnames <- c()

for (i in temp_colnames) {
  if (str_detect(i, "-11A-") == TRUE) {
    print("OK") 
    temp_column_name = paste0(i, "_colon")
    new_temp_colnames = c(new_temp_colnames, temp_column_name)
  } else if ((str_detect(i, "-01A-") == TRUE)) {
    print("Also OK") 
    temp_column_name = paste0(i, "_cancer")
    new_temp_colnames = c(new_temp_colnames, temp_column_name)
  } else {
    print("NOT OK")
    new_temp_colnames = c(new_temp_colnames, i)
  } 
}

colnames(patients_lateStages) <- new_temp_colnames

patients_late <- inner_join(patients_lateStages, patients_GSE50760, by = "Symbol")


### MODIFY COLUMNS OF COUNT MATRIX FOR EARLY STAGES

patients_earlyStages <- readRDS("~/Desktop/HT projects/CRC_transcriptomics_proj/data/TCGA/Illumina_HiSeq_Results/patients_Earlystages.RDS") %>%
  as.data.frame()
dim(patients_earlyStages)
patients_earlyStages <- patients_earlyStages[early_down$Symbol,]
patients_earlyStages <- rownames_to_column(patients_earlyStages, "Symbol")

temp_colnames <- colnames(patients_earlyStages)
new_temp_colnames <- c()

for (i in temp_colnames) {
  if (str_detect(i, "-11A-") == TRUE) {
    print("OK") 
    temp_column_name = paste0(i, "_colon")
    new_temp_colnames = c(new_temp_colnames, temp_column_name)
  } else if ((str_detect(i, "-01A-") == TRUE)) {
    print("Also OK") 
    temp_column_name = paste0(i, "_cancer")
    new_temp_colnames = c(new_temp_colnames, temp_column_name)
  } else {
    print("NOT OK")
    new_temp_colnames = c(new_temp_colnames, i)
  } 
}

colnames(patients_earlyStages) <- new_temp_colnames

patients_early <- patients_earlyStages

rm(patients_GSE50760, patients_lateStages, patients_earlyStages)

##### Normalize + z-score

### Early stages
p <- patients_early
HeatmapData <- dplyr::mutate_at(p, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% DEA_Early_Filtered$Symbol)
HeatmapData_input <- column_to_rownames(HeatmapData_input, "Symbol")
HeatmapData_input <- as.data.frame(t(HeatmapData_input))
HeatmapData_input <- rownames_to_column(HeatmapData_input, var = "barcodes")
HeatmapData_input$Project <- with(HeatmapData_input,
                                  ifelse(test = grepl("TCGA", HeatmapData_input$barcodes),
                                         yes = "TCGA", 
                                         no = "GSE50760"))
HeatmapData_input$Status <- with(HeatmapData_input,
                                 ifelse(test = grepl("_cancer", HeatmapData_input$barcodes),
                                         yes = "Cancer", 
                                        no = "Colon"))
HeatmapData_input <- HeatmapData_input[, c(ncol(HeatmapData_input),
                                           (ncol(HeatmapData_input)-1),
                                           1:(ncol(HeatmapData_input)-2))]
HeatmapData_matrix <- as.matrix(HeatmapData_input[, c(4:ncol(HeatmapData_input))])
dim(HeatmapData_matrix)

sample.colors = c("gold", "black")
sample.colors.2 = c("dodgerblue", "deeppink2")
names(sample.colors) = c("Colon", "Cancer")
names(sample.colors.2) = c("TCGA", "GSE50760")
status_info = data.frame(Status = HeatmapData_input$Status)
status_info.2 = data.frame(Project = HeatmapData_input$Project)

fontsize = 0.6

stage_early <- readRDS(file.path(data.dir, "unique_4_early.RDS")) # stage_early.RDS unique_4_early.RDS, stage_early_down.RDS


flog.debug("Print column names")
ha_number <- which(colnames(HeatmapData_matrix) %in% stage_early) # stage_early$Symbol
length(ha_number)

flog.debug("Print heatmap annotation")
ha = rowAnnotation(foo = anno_mark(at = ha_number, labels = stage_early)) # stage_early$Symbol

heatmap_final <-
  Heatmap(t(HeatmapData_matrix),
          cluster_columns = TRUE,
          column_names_side = "top",
          column_dend_side = "top",
          column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
          row_names_side = "left",
          #row_dend_side = "left",
          show_row_names = FALSE,
          row_names_gp = gpar(cex = fontsize),
          row_dend_width = unit(3, "cm"),
          clustering_distance_rows = "euclidean",
          clustering_method_rows = "ward.D",
          name = "z-score",
          # column_title = "Early stages of cancerogenesis - Cancer vs. Colon",
          column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
          right_annotation = ha,
          row_title_rot = 0,
          col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
          row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
          top_annotation = HeatmapAnnotation(
            Status = status_info$Status,
            Project = status_info.2$Project,
            col = list(Status = c("Cancer" = "gold", "Colon" = "black"),
                       Project = c("TCGA" = "dodgerblue", "GSE50760" = "deeppink2")),
            show_legend = TRUE,
            show_annotation_name = FALSE,
            annotation_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                                           legend_direction = "horizontal",
                                           title_position = "topcenter")
            # annotation_name_side = "left"
          ),
          heatmap_legend_param = list(
            title = "z-score", 
            title_gp = gpar(fontsize = 14, fontface = "bold"),
            title_position = "leftcenter-rot",
            direction = "vertical"
          ))

draw(heatmap_final, heatmap_legend_side = "right", 
     annotation_legend_side = "top", merge_legend = TRUE)


### Late stages

p <- patients_late
p <- dplyr::filter(p, p$Symbol %in% DEA_Late_Filtered$Symbol)

# TCGA
p_TCGA <-  p[, 1:25]
HeatmapData_TCGA <- dplyr::mutate_at(p_TCGA, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData_TCGA[, c(2:ncol(HeatmapData_TCGA))], 1, mean)
StDev <- apply(HeatmapData_TCGA[, c(2:ncol(HeatmapData_TCGA))], 1, sd)
HeatmapData_TCGA_input <- mutate_at(HeatmapData_TCGA, vars(-c(Symbol)), funs(((.)-Av)/StDev))

# GSE50760
p_GSE <-  p[, c(1,26:61)]
HeatmapData_GSE <- dplyr::mutate_at(p_GSE, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData_GSE[, c(2:ncol(HeatmapData_GSE))], 1, mean)
StDev <- apply(HeatmapData_GSE[, c(2:ncol(HeatmapData_GSE))], 1, sd)
HeatmapData_GSE_input <- mutate_at(HeatmapData_GSE, vars(-c(Symbol)), funs(((.)-Av)/StDev))

# bind
HeatmapData_input <- cbind(HeatmapData_TCGA_input, HeatmapData_GSE_input[, c(2:ncol(HeatmapData_GSE_input))])

HeatmapData_input <- column_to_rownames(HeatmapData_input, "Symbol")
HeatmapData_input <- as.data.frame(t(HeatmapData_input))
HeatmapData_input <- rownames_to_column(HeatmapData_input, var = "barcodes")
HeatmapData_input$Project <- with(HeatmapData_input,
                                  ifelse(test = grepl("TCGA", HeatmapData_input$barcodes),
                                         yes = "TCGA", 
                                         no = "GSE50760"))
HeatmapData_input$Status <- with(HeatmapData_input,
                                 ifelse(test = grepl("_cancer", HeatmapData_input$barcodes),
                                        yes = "Cancer", 
                                        no = "Colon"))
HeatmapData_input <- HeatmapData_input[, c(ncol(HeatmapData_input),
                                           (ncol(HeatmapData_input)-1),
                                           1:(ncol(HeatmapData_input)-2))]
HeatmapData_matrix <- as.matrix(HeatmapData_input[, c(4:ncol(HeatmapData_input))])

sample.colors = c("gold", "black")
sample.colors.2 = c("dodgerblue", "deeppink2")
names(sample.colors) = c("Cancer", "Colon")
names(sample.colors.2) = c("TCGA", "GSE50760")
status_info = data.frame(Status = HeatmapData_input$Status)
status_info.2 = data.frame(Project = HeatmapData_input$Project)

fontsize = 0.6

stage_late <- readRDS(file.path(data.dir, "unique_4_late.RDS")) # stage_late.RDS, unique_4_late.RDS, stage_late_down.RDS

flog.debug("Print column names")
ha_number <- which(colnames(HeatmapData_matrix) %in% stage_late) # stage_late$Symbol
flog.debug("Print heatmap annotation")
ha = rowAnnotation(foo = anno_mark(at = ha_number, labels = stage_late)) # stage_late$Symbol


heatmap_final <-
  Heatmap(t(HeatmapData_matrix),
          cluster_columns = TRUE,
          column_names_side = "top",
          column_dend_side = "top",
          column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
          row_names_side = "left",
          #row_dend_side = "left",
          show_row_names = FALSE,
          row_names_gp = gpar(cex = fontsize),
          row_dend_width = unit(3, "cm"),
          clustering_distance_rows = "euclidean",
          clustering_method_rows = "ward.D",
          name = "z-score",
          # column_title = "Early stages of cancerogenesis - Cancer vs. Colon",
          column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
          right_annotation = ha,
          row_title_rot = 0,
          col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
          row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
          top_annotation = HeatmapAnnotation(
            Status = status_info$Status,
            Project = status_info.2$Project,
            col = list(Status = c("Cancer" = "gold", "Colon" = "black"),
                       Project = c("TCGA" = "dodgerblue", "GSE50760" = "deeppink2")),
            show_legend = TRUE,
            show_annotation_name = FALSE,
            annotation_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                                           legend_direction = "horizontal",
                                           title_position = "topcenter")
            # annotation_name_side = "left"
          ),
          heatmap_legend_param = list(
            title = "z-score", 
            title_gp = gpar(fontsize = 14, fontface = "bold"),
            title_position = "leftcenter-rot",
            direction = "vertical"
          ))

draw(heatmap_final, heatmap_legend_side = "right", 
     annotation_legend_side = "top", merge_legend = TRUE)



flog.debug("Boxplot")

GOI <- c("TSG101", "VPS28", "VPS37A", "VPS37B", "VPS37C", "VPS37D")

TPM_calculator(patients_early)

ESCRTs_earlyStage <- filter(TPM, TPM$Symbol %in% GOI)
ESCRTs_earlyStage <- column_to_rownames(ESCRTs_earlyStage, "Symbol")
ESCRTs_earlyStage <- t(ESCRTs_earlyStage) %>% as.data.frame()
ESCRTs_earlyStage <- rownames_to_column(ESCRTs_earlyStage, "sampleID")
ESCRTs_earlyStage$Status <- factor(ifelse(
  test = str_detect(ESCRTs_earlyStage$sampleID, "cancer"),
  yes = "Cancer",
  no = "Colon"
), levels = c("Colon", "Cancer"))

bp_VPS37A_early = ggplot(
  data = ESCRTs_earlyStage,
  mapping = aes(x = Status, y = ESCRTs_earlyStage$VPS37A, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37A\n Early stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37A_early)

bp_VPS37B_early = ggplot(
  data = ESCRTs_earlyStage,
  mapping = aes(x = Status, y = ESCRTs_earlyStage$VPS37B, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37B\n Early stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37B_early)

bp_VPS37C_early = ggplot(
  data = ESCRTs_earlyStage,
  mapping = aes(x = Status, y = ESCRTs_earlyStage$VPS37C, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37C\n Early stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37C_early)

TPM_calculator(patients_late)

ESCRTs_lateStage <- filter(TPM, TPM$Symbol %in% GOI)
ESCRTs_lateStage <- column_to_rownames(ESCRTs_lateStage, "Symbol")
ESCRTs_lateStage <- t(ESCRTs_lateStage) %>% as.data.frame()
ESCRTs_lateStage <- rownames_to_column(ESCRTs_lateStage, "sampleID")
ESCRTs_lateStage$Status <- factor(ifelse(
  test = str_detect(ESCRTs_lateStage$sampleID, "cancer"),
  yes = "Cancer",
  no = "Colon"
), levels = c("Colon", "Cancer"))

bp_VPS37A_late = ggplot(
  data = ESCRTs_lateStage,
  mapping = aes(x = Status, y = ESCRTs_lateStage$VPS37A, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37A\n Advanced stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37A_late)

bp_VPS37B_late = ggplot(
  data = ESCRTs_lateStage,
  mapping = aes(x = Status, y = ESCRTs_lateStage$VPS37B, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37B\n Advanced stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37B_late)

bp_VPS37C_late = ggplot(
  data = ESCRTs_lateStage,
  mapping = aes(x = Status, y = ESCRTs_lateStage$VPS37C, fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  # ggtitle("Expression of VPS37C\n Advanced stages of carcinogenesis") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(bp_VPS37C_late)