rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(tidyverse)
  library(cowplot)
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
out.dir = file.path(proj.dir, "data", "output", "TCGA", "RNASeq")
dir.create(out.dir, recursive = TRUE)


flog.debug("Load files")

Viz_Earlystages <- readRDS(file.path(data.dir, "Viz_input_Earlystages.RDS"))
Viz_Latestages <- readRDS(file.path(data.dir, "Viz_input_Latestages.RDS"))
Stats_Earlystages <- readRDS(file.path(data.dir, "Stat_Analysis_Earlystages.RDS"))
Stats_Latestages <- readRDS(file.path(data.dir, "Stat_Analysis_Latestages.RDS"))


flog.debug("Prepare boxplots")

GOI <- c("VPS37A", "VPS37B", "VPS37C")

viz_input = Viz_Earlystages
# viz_input = Viz_Latestages
viz_input = viz_input[GOI]

list.plots <- NULL

for (i in 1:length(viz_input)) {
  gene_name = names(viz_input[i])
  temp_data = viz_input[[i]]
  bp = ggplot(
    data = temp_data,
    mapping = aes(x = Status, y = temp_data[,gene_name], fill = Status)) +
    geom_boxplot() +
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
          panel.grid.minor = element_blank())
  print(bp)
  list.plots[[gene_name]] = bp
  assign(paste0("Boxplot_", gene_name), bp)
}


# get legend
bp <- ggplot(
  data = temp_data,
  mapping = aes(x = Status, y = temp_data[,gene_name], fill = Status)) +
  geom_boxplot() +
  geom_dotplot(colour = "black", fill = "black", binaxis='y', stackdir = 'center', dotsize = 1) +
  scale_fill_manual(values = c("dodgerblue", "deeppink2")) +
  ggtitle(paste0(gene_name)) +
  theme(plot.title = element_text(face = "bold")) +
  xlab("Disease status") + ylab("Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 21, colour = 'black', face = "bold"),
        axis.title.x = element_text(face = "bold", colour = 'black', size = 14), 
        axis.title.y = element_text(face = "bold", colour = 'black', size = 14), 
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
legend <- get_legend(bp)

title = ggdraw() + 
  draw_label(
    paste0("TCGA"),
    fontface = 'bold',
    size = 21,
    x = 0.5,
    hjust = 0.5
  ) +
  theme()

rm(bp)

# fn = "GOI_earlyStages"
fn = "GOI_lateStages"
png(filename = file.path(out.dir, paste0(format(Sys.time(),  "%Y-%m-%d"), fn, ".png")),
    width = 1400,
    height = 800,
    units = "px")
plot <- plot_grid(
  title,
  plot_grid(plotlist = list.plots,
            nrow = 2, ncol = 2,
            rel_widths = c(1, 1),
            rel_heights = c(1, 1)),
  legend,
  nrow = 3,
  rel_heights = c(0.1, 1, 0.1)) 
print(plot)
dev.off()


sessionInfo()