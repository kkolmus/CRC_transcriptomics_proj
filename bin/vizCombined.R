rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scran")

suppressPackageStartupMessages({
  library(futile.logger)
  library(tidyverse)
  library(plyr)
  library(SingleCellExperiment)
  library(scran)
  library(ggforce)
})

flog.threshold(DEBUG)


flog.debug("Set project environment")

proj.dir = "~/Desktop/HT projects/CRC_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()
TCGA.data.dir = file.path(proj.dir, "data", "TCGA", "RNASeq")
GSE.data.dir = file.path(proj.dir, "data", "input", "GSE50760")


flog.debug("Perform statistical analysis on individual datasets")

# Genes of interest

GOI <- c("VPS37A", "VPS37B", "VPS37C")

# TCGA

Viz_Earlystages <- readRDS(file.path(TCGA.data.dir, "Viz_input_Earlystages.RDS"))
Viz_Earlystages <- Viz_Earlystages[GOI]
Viz_Earlystages <- ldply(Viz_Earlystages, data.frame)
Viz_Earlystages <- Viz_Earlystages[, c(1,2,4,3,5,6)]
 
Viz_Latestages <- readRDS(file.path(TCGA.data.dir, "Viz_input_Latestages.RDS"))
Viz_Latestages <- Viz_Latestages[GOI]
Viz_Latestages <- ldply(Viz_Latestages, data.frame)
Viz_Latestages <- Viz_Latestages[, c(1,2,4,3,5,6)]

for (i in 1:length(GOI)) {
  gene = GOI[i]
  col_number = 3+i
  TCGA_CRC = filter(Viz_Latestages, 
                     Viz_Latestages$`.id` == gene &
                       Viz_Latestages$Status == "Cancer")
  TCGA_CRC_mean = mean(TCGA_CRC[, col_number])
  assign(paste0(gene, "_TCGA_CRC"), TCGA_CRC_mean)
  
  TCGA_Colon = filter(Viz_Latestages, 
                        Viz_Latestages$`.id` == gene &
                          Viz_Latestages$Status == "Colon")
  TCGA_Colon_mean = mean(TCGA_Colon[, col_number])
  assign(paste0(gene, "_TCGA_Colon"), TCGA_Colon_mean)
  
  TCGA_stats = t.test(
    x = TCGA_CRC[, col_number],
    y = TCGA_Colon[, col_number],
    alternative = "two.sided",
    paired = TRUE
  )
  TCGA_stats = TCGA_stats[["p.value"]]
  assign(paste0(gene, "_TCGA_pval"), TCGA_stats)
  
  rm(TCGA_CRC, TCGA_Colon, TCGA_CRC_mean, TCGA_Colon_mean, TCGA_stats)
}


# GSE50670

Viz_GSE50760 <- readRDS(file.path(GSE.data.dir, "Viz_input_GSE50760.RDS"))
Viz_GSE50760 <- Viz_GSE50760[GOI]
Viz_GSE50760 <- ldply(Viz_GSE50760, data.frame)
Viz_GSE50760 <- Viz_GSE50760[, c(1,2,4,3,5,6)]
Viz_GSE50760$Status <- ifelse(
  test = Viz_GSE50760$Status %in% "Colorectal cancer",
  yes = "Cancer",
  no = "Colon"
)

for (i in 1:length(GOI)) {
  gene = GOI[i]
  col_number = 3+i
  GSE50760_CRC = filter(Viz_GSE50760, 
                        Viz_GSE50760$`.id` == gene &
                      Viz_GSE50760$Status == "Cancer")
  GSE50760_CRC_mean = mean(GSE50760_CRC[, col_number])
  assign(paste0(gene, "_GSE50760_CRC"), GSE50760_CRC_mean)
  
  GSE50760_Colon = filter(Viz_GSE50760, 
                          Viz_GSE50760$`.id` == gene &
                            Viz_GSE50760$Status == "Colon")
  GSE50760_Colon_mean = mean(GSE50760_Colon[, col_number])
  assign(paste0(gene, "_GSE50760_Colon"), GSE50760_Colon_mean)
  
  GSE50760_stats = t.test(
    x = GSE50760_CRC[, col_number],
    y = GSE50760_Colon[, col_number],
    alternative = "two.sided",
    paired = TRUE
  )
  GSE50760_stats = GSE50760_stats[["p.value"]]
  assign(paste0(gene, "_GSE50760_pval"), GSE50760_stats)
  
  rm(GSE50760_CRC, GSE50760_Colon, GSE50760_CRC_mean, GSE50760_Colon_mean, GSE50760_stats)
}


flog.debug("Perform meta-analysis")

# VPS37A
Av_VPS37A_Colon <- mean(c(VPS37A_GSE50760_Colon, VPS37A_TCGA_Colon))
Av_VPS37A_CRC <- mean(c(VPS37A_GSE50760_CRC, VPS37A_TCGA_CRC))

StDev_VPS37A_Colon <- combinePValues(VPS37A_GSE50760_pval, VPS37A_TCGA_pval,
                                     method = "fisher")

# VPS37B
Av_VPS37B_Colon <- mean(c(VPS37B_GSE50760_Colon, VPS37B_TCGA_Colon))
Av_VPS37B_CRC <- mean(c(VPS37B_GSE50760_CRC, VPS37B_TCGA_CRC))

StDev_VPS37B_Colon <- combinePValues(VPS37B_GSE50760_pval, VPS37B_TCGA_pval,
                                     method = "fisher")

# VPS37C
Av_VPS37C_Colon <- mean(c(VPS37C_GSE50760_Colon, VPS37C_TCGA_Colon))
Av_VPS37C_CRC <- mean(c(VPS37C_GSE50760_CRC, VPS37C_TCGA_CRC))

StDev_VPS37C_Colon <- combinePValues(VPS37C_GSE50760_pval, VPS37C_TCGA_pval,
                                     method = "fisher")


flog.debug("Data vizualization")

lateStages <- rbind(Viz_GSE50760, Viz_Latestages)
lateStages$Status <- factor(lateStages$Status, level = c("Colon", "Cancer"))

for (i in 1:length(GOI)) {
  gene_name = GOI[i]
  temp_data = filter(lateStages, lateStages$`.id` == gene_name)
  bp = ggplot(
    data = temp_data,
    mapping = aes(x = Status, y = temp_data[,gene_name], fill = Status)) +
    geom_violin() +
    geom_dotplot(colour = "black", fill = "black", binaxis = 'y', stackdir = 'center', dotsize = 1) +
    scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
    ggtitle(paste0(gene_name)) +
    theme(plot.title = element_text(face = "bold")) +
    xlab("Disease status") + ylab("Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
          axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
          axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
          axis.text = element_text(size = 14, color = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    facet_zoom(
      ylim = c(0, 100),
      zoom.size = 1)
  print(bp)
  list.plots[[gene_name]] = bp
  assign(paste0("ViolinPlot_", gene_name), bp)
}

sessionInfo()