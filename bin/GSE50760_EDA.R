### LIBS ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(tidyverse)
  library(cowplot)
})

### MAIN ###

flog.threshold(NULL)

flog.debug("Prepare working environment")

rm(list = ls())

projID =  "GSE50760"

proj.dir <- "~/Desktop/HT projects/CRC_transcriptomics_proj"
data.dir <- file.path(proj.dir, "data")
input.dir <- file.path(data.dir, "input", projID)
output.dir <- file.path(data.dir, "output", projID)
dir.create(output.dir, recursive = TRUE)


flog.debug("Load data")

pdata <- readRDS(file.path(input.dir, "pdata.RDS"))
pdata$status <- c(rep("CRC", 18), rep("Colon", 18))
pdata$newIDs <- paste(pdata$status, pdata$patient_id, sep = "_") 
 
n.counts <- readRDS(file.path(input.dir, "n.counts.RDS"))
colnames(n.counts) <- c("Symbol", pdata$newIDs)

genes <- read.table(file.path(data.dir, "endocytic_genes.txt")) %>% t() %>% c()
genes <- substring(genes, 1, (nchar(genes)-1))
length(genes)
# 445

flog.debug("Select endocytic genes")

n.counts_subset <- filter(
  n.counts, 
  n.counts$Symbol %in% genes
)

genes_subset <- n.counts_subset$Symbol
length(genes_subset)
# 435


flog.debug("Convert df into a list with separate slots for each gene")

viz_input <- setNames(vector("list", length(genes_subset)), genes_subset)

for(i in 1:nrow(n.counts_subset)) {
  print(n.counts_subset[i,1])
  temp_gene = n.counts_subset[i,1]
  temp_data = n.counts_subset[i,]
  temp_data = as.data.frame(t(temp_data))
  temp_data = rownames_to_column(temp_data, "metadata")
  temp_column = c("newID", temp_gene)  # as.character(unname(unlist(as.data.frame(t(temp_data[1,])))))
  colnames(temp_data) = temp_column
  temp_data = temp_data[-1,]
  rownames(temp_data) = NULL
  temp_data$Status <- factor(ifelse(
    test = str_detect(temp_data$newID, "CRC"),
    yes = "Colorectal cancer",
    no = "Colon"),
    levels = c("Colon", "Colorectal cancer"))
  temp_data[,2] = as.numeric(as.character(temp_data[,2]))
  viz_input[[temp_gene]] = temp_data
}

rm(temp_column, temp_gene, temp_data, i)

saveRDS(viz_input, file.path(input.dir, paste0("Viz_input_", projID, ".RDS")))


flog.debug("Statsitical analysis")

stats <- setNames(vector("list", length(genes_subset)), genes_subset)

for (i in 1:length(viz_input)) {
  print(names(viz_input[i]))
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

rm(stats, test, temp_data, i)

saveRDS(stats_df_pval_subset, file.path(input.dir, paste0("Stats", projID, ".RDS")))


flog.debug("Prepare boxplots")

GOI <- c("VPS37A", "VPS37B", "VPS37C")

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

rm(bp, i)

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
    paste0(projID),
    fontface = 'bold',
    size = 21,
    x = 0.5,
    hjust = 0.5
  ) +
  theme()

rm(bp)

fn = "GOI"
png(filename = file.path(output.dir, paste0(format(Sys.time(),  "%Y-%m-%d"), fn, ".png")),
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