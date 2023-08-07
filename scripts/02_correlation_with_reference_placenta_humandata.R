# This script outputs plots of the correlation with reference placenta clusters: paper Fig 1C
# code adapted from https://github.com/archavan/covid-placenta

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  
})

source('scripts/00_get_started.R')

#  load data ----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE.rds"
so = readRDS(so_file)
Idents(so) = "sub_cluster_names"
DefaultAssay(so) <- "RNA"

# Reference datasets ----
ref_file="/cluster/tufts/slonimlab/rbator01/reference_data/three_dataset_merged_markers/vento-surya-covid_control-joined.csv"
ref <- read.csv(ref_file, stringsAsFactors = FALSE)

# rename
names(ref) <- gsub("lu", "l", names(ref))
names(ref) <- gsub("surya", "s", names(ref))
names(ref) <- gsub("vento", "v", names(ref))

# Correlation matrix
so.avg <- AverageExpression(so)
so.avg.rna <- so.avg$RNA
so.avg.rna = data.frame(so.avg.rna)
colnames(so.avg.rna) = gsub("^","P_",colnames(so.avg.rna))
so.avg.rna$gene_name <- rownames(so.avg.rna)

## join reference and placenta
plac_ref=dplyr::inner_join(so.avg.rna, ref, by = "gene_name")
plac_ref$HGNC.symbol = NULL

# Annotation by correlation ----

refdata <- c("l","s","v")

# build correlation matrix from expression data
cor.matrix <- cor(plac_ref[, names(plac_ref)[names(plac_ref) != "gene_name" ]], 
                  method = 'spearman')

# reorder correlation matrix based on clustering
dd <- as.dist((1 - cor.matrix))
hc <- hclust(dd, method = 'complete')
cor.matrix <- cor.matrix[hc$order, hc$order]

# melt correlation matrix
cormat <- reshape2::melt(cor.matrix, na.rm = T)
cormat$Var1_source <- sapply(strsplit(as.character(cormat$Var1), split = "_"), "[[", 1)
cormat$Var2_source <- sapply(strsplit(as.character(cormat$Var2), split = "_"), "[[", 1)

# subset
cormat <- cormat[which(cormat$Var1_source %in% refdata & 
                         cormat$Var2_source %in% c("P")), ]

# top 3 matches with highest correlation coefficients
topmatch <- data.frame( # empty df to fill top matches from the loop below
  cluster = unique(as.character(cormat$Var2)) %>% sort(),
  vento.top1 = NA,
  surya.top1 = NA,
  covid.top1 = NA,
  vento.top2 = NA,
  surya.top2 = NA,
  covid.top2 = NA,
  vento.top3 = NA,
  surya.top3 = NA,
  covid.top3 = NA
)

cormat$top3 <- NA

for(i in unique(cormat$Var2)){
  for(j in refdata){
    # identify top3 match indices
    ind <- which(cormat$Var2 == i & cormat$Var1_source == j)
    val <- cormat$value[ind]
    top3ind <- ind[order(val, decreasing = TRUE)[1:3]]
    
    # assign match ranking to correlation data for plotting
    cormat$top3[top3ind[1]] <- "1"
    cormat$top3[top3ind[2]] <- "2"
    cormat$top3[top3ind[3]] <- "3"
    
    # assign top matches to topmatches dataframe
    topmatch[topmatch$cluster == i, paste0(j, ".top1")] <- gsub(paste0(j, "_"), "", as.character(cormat$Var1[top3ind[1]]))
    topmatch[topmatch$cluster == i, paste0(j, ".top2")] <- gsub(paste0(j, "_"), "", as.character(cormat$Var1[top3ind[2]]))
    topmatch[topmatch$cluster == i, paste0(j, ".top3")] <- gsub(paste0(j, "_"), "", as.character(cormat$Var1[top3ind[3]]))
  }
}
library(viridis)

# Plot all, this is not the most informative plot, it is easier to see if you break down the individual datasets as below -----

# output directory for plots
p <- ggplot(data = cormat,
            aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_viridis(name = "Spearman\nCorrelation") +
  geom_point(aes(Var1, Var2, alpha = top3),
             size = 1.5, shape = 19, stroke  = 0) +
  scale_alpha_manual(values = c(1, 0.5, 0.25),
                     breaks = c(1, 2, 3),
                     name = "top3\nwithin dataset", na.value = 0) +
  coord_fixed(ratio = 1) +
  xlab("reference") +
  ylab("query") +
  labs(caption = "For each celltype in query, black points represent top 3 celltypes from reference with highest correlation.",
       title = paste0("pl vs.reference")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8,
                                   hjust=1, vjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.ticks.length = unit(0.15, units = c('lines')),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        plot.caption = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.8, units = c('lines')))

show(p)
ggsave(p, filename = "analysis/plots/ref_compare/pl_subset_clusters_res_0.2_vs_reference_corrplot_clusters.pdf",
       device = "pdf", width = 15, height = 7, units = "in")


# plotting function
clustAnnoPlot <- function(dat, query, reference, plot_title) {
  
  dat <- dat[which(dat$Var2_source %in% query & dat$Var1_source == reference), ]

  p <- ggplot(data = dat, 
              aes(Var1, Var2, fill = value)) +
    geom_tile(colour = "white") +
    scale_fill_viridis(name = "Spearman\nCorrelation") +
    geom_point(aes(Var1, Var2, alpha = top3),
               size = 1.5, shape = 19, stroke  = 0) +
    scale_alpha_manual(values = c(1, 0.5, 0.25), 
                       breaks = c(1, 2, 3),
                       name = "top3", na.value = 0) +
    coord_fixed(ratio = 1) +
    xlab("reference") +
    ylab("query") +
    
    labs(title = plot_title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=1, vjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.ticks.length = unit(0.15, units = c('lines')),
          legend.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          plot.caption = element_text(size = 7),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.key.size = unit(0.8, units = c('lines')))
  
  return(p)
}

clustAnnoPlot_notop <- function(dat, query, reference, plot_title) {
  
  dat <- dat[which(dat$Var2_source %in% query & dat$Var1_source == reference), ]
  
  p <- ggplot(data = dat, 
              aes(Var1, Var2, fill = value)) +
    geom_tile(colour = "white") +
    scale_fill_viridis(name = "Spearman\nCorrelation") +
    coord_fixed(ratio = 1) +
    xlab("reference") +
    ylab("query") +
    labs(title = plot_title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=1, vjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.ticks.length = unit(0.15, units = c('lines')),
          legend.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          plot.caption = element_text(size = 7),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.key.size = unit(0.8, units = c('lines')))
  
  return(p)
}

# Vento-Tormo et al
anno.vento <- clustAnnoPlot(dat = cormat, query = "P", reference = "v", plot_title="Vento-Tormo et al")
cowplot::ggsave2(anno.vento, device = "pdf", width = 8, height = 4, units = "in",
                 filename = paste0(outdir, "vento_subset_heatmap.pdf"))

#"Lu-Culligan et al
anno.covid <- clustAnnoPlot(dat = cormat, query = c("P"), reference = "l", plot_title="Lu-Culligan et al")
cowplot::ggsave2(anno.covid, device = "pdf", width = 6, height = 4, units = "in",
                 filename = "plots_compare_reference/lu_c_subset_heatmap.pdf")
show(anno.covid)
#Suryavanshi et al
anno.surya <- clustAnnoPlot(dat = cormat, query = c("P"), reference = "s", plot_title="Suryawanshi et al")
cowplot::ggsave2(anno.surya, device = "pdf", width = 6, height = 4, units = "in",
                 filename = "plots_compare_reference/surya_subset_heatmap.pdf")
show(anno.surya)
