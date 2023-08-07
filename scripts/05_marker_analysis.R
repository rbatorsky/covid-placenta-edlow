# Plot markers, deg and enrichment results: Paper Fig 1B,D,E

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(openxlsx)
  library(here)
})

source('scripts/00_get_started.R')

# load data   ----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE_norm.rds"
so = readRDS(so_file)
res="sub_cluster_names"
Idents(so) = res
DefaultAssay(so) <- "RNA"

# sex marker genes ----
Idents(so) = "covid_status"
DimPlot(so, split.by="covid_status")

# Sex markers: Fig 1B----
sex_genes=rev(c("KDM5D","UTY","TSIX","XIST","EIF2S3","DDX3Y"))
sex_genes_few=rev(c("XIST","FLNA","PAGE4","DDX3Y","GPR34","RPS4Y1"))
sex_genes_few=rev(c("XIST","DDX3Y"))

so_m = subset(so, sex == "M")
#so_m = subset(so_m, donor != "D145")
DefaultAssay(so_m) <- "RNA"

p0 = DotPlot(so_m, features=sex_genes_few, 
        scale = F) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_fill_binned(type = "viridis")

ggsave(p0, filename = "analysis/markers/plots/sex_markers_dot_hbc_subset_clusters.pdf", 
       useDingbats = F, height = 3, width = 6, units = "in")

dev.off()

# Find some genes on ChrX ChrY and average, not used in paper ----
gtf = read.table("/cluster/tufts/bio/data/genomes/HomoSapiens/UCSC/hg38/Annotation/Genes/genes.gtf", header = FALSE, sep = '\t') 
chrxy = gtf %>% dplyr::filter(V1 %in% c("chrX", "chrY"))
sep = separate(data = chrxy, col = V9, into = c("A","B","C","D"), sep = ";")
sep = sep %>% dplyr::select(c("V1","A"))
sep$A = gsub("gene_id ","", sep$A)
xygenes = unique(sep$A)

# Average the genes
so.avg <- AverageExpression(so)
counts = data.frame(so.avg$RNA)
counts$gene = rownames(counts)
counts_xy = counts %>% dplyr::filter(gene %in% xygenes)
nrow(counts_xy)
view(counts_xy)

counts_xy_select = counts_xy %>% dplyr::select(c("HBCs_1","PAMs","gene"))

counts_xy_select_filter = counts_xy_select %>% dplyr::filter(HBCs_1 <1 & PAMs > 1 | PAMs<1 & HBCs_1>1)

ggplot(counts_xy_select_filter, aes(x = HBCs_1, y = PAMs, label = gene)) + 
  geom_point() + 
  xlim(c(0,5)) + 
  ylim(c(0,5)) + 
  geom_text()
  

# Read the markers ----
markers = read.csv("analysis/markers/markers_seuratobj_subset_integrated_HBCs_USETHISONE.csv")
markers %>% dplyr::filter(cluster == "HBCs_2")

markers_filter = markers %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.25)
table(markers_filter$cluster)


top10.sub <- markers_filter %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)


top5.sub <- markers_filter %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  arrange(factor(cluster, levels=c("HBCs_0","HBCs_1","HBCs_2","HBCs_3","HBCs_4","HBCs_5","HBCs_6","HBCs_7","PAMs","Monocytes")))

so.avg <- AverageExpression(so, return.seurat = TRUE)

#  averaged heatmap -----
p1 = DoHeatmap(so.avg,features=top5.sub$gene, raster = FALSE) + 
  theme(text = element_text(size = 12)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256), na.value = "white")

show(p1)
dev.off()
ggsave(p1, filename = "analysis/markers/plots/marker_top5_heatmap_hbc_subset_clusters.pdf", 
       useDingbats = F, height = 8, width = 7, units = "in")

# read the GO enrichment results ----
ck = readRDS("analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25.rds")

select_terms_workbook = read.xlsx("analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25.xlsx")

select_terms_workbook = select_terms_workbook %>%
  dplyr::filter(old_select == 1) %>%
  mutate(category = factor(category, levels = c("immune", "metabolic", "stress", "vascular", "ribo", "protein folding", "migration", "phago", "actin"))) %>%
  arrange(category)

marker_term_select_combine = unique(select_terms_workbook$Description)

# Make an ordered GO dot plot: Fig 1E ----

p = ordered_cp_dotplot(ck, marker_term_select_combine)

print(p)
ggsave(p, filename = "markers/plots/select_terms_go.pdf",
       device = cairo_pdf,
       width = 9, height = 12, 
       units = "in")


# non-ordered dot plot, not used in paper ----
p = clusterProfiler::dotplot(ck_select, x =~cluster,showCategory =50, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(p, filename = "deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_top20.pdf",
       device = cairo_pdf,
       width = 13, height = 18,
       units = "in")

# bar chart visual: not used ----

ck_top  = data.frame(ck)  %>% 
  arrange(p.adjust) %>% group_by(cluster) %>% 
  top_n(10, wt=-p.adjust)

ck_top$nlogp = -log10(ck_top$p.adjust)
unique(ck_top$cluster)

color=sample(brewer.pal(8, "Spectral"))
categories=unique(ck_top$Cluster)
names(color) = categories
outdir="/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/markers/"

p0 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_0"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47) + xlab("")

p1 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_1"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47)+ ylab("")+ xlab("")


p2 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_2"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47)+ ylab("")+ xlab("")

p3 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_3"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47) + ylab("")+ xlab("")

p5 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_5"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47) + ylab("")+ xlab("")

ppam = ggplot(ck_top %>% dplyr::filter(cluster=="PAMs"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic() + 
  ylim(0,47)+ ylab("")+ xlab("")

p7 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_7"),
              aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,47)+ ylab("")+ xlab("")

ggarrange(p7 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 105) ),
          ppam +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l=0) ), 
          p5 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 0) ), 
          p3 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 90)),
          p2 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 105) ),
          p1 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 105)),
          p0 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 105) ),
          ncol = 1)


for (i in unique(ck_top$cluster)){
  
  ck_top_1 = ck_top %>% dplyr::filter(Cluster %in% c(i))
  p = ggbarplot(ck_top_1, x = "Description", y = "nlogp",
                fill = "Cluster",
                sort.val = "asc",  
                sort.by.groups = FALSE,
                ylab = "-log(p.adj)") + 
    theme(text = element_text(size = 12)) + coord_flip() + 
    theme(legend.position="right") + 
    facet_grid(~cluster, scales = "free_y") + 
    scale_y_continuous(expand = c(0,0)) 
  
  ggsave(p, filename = paste0(outdir,paste0(i,"_gobar.pdf")),
         device = cairo_pdf,
         width = 9, height = 3,
         units = "in")
  
}

# all clusters, before subset - not used ----

so_all_markers = read.csv("/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/markers/so_all_markers_res_0.3.csv")
head(so_all_markers)

top10.sub <- so_all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

DotPlot(so_all, features=unique(top10.sub$gene)) +
 theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8)) + 
  coord_flip()

so_all.avg <- AverageExpression(so_all, return.seurat = TRUE)


p1 = DoHeatmap(so_all.avg,features=top10.sub$gene, raster = FALSE) + 
  theme(text = element_text(size = 10)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))

show(p1)
ggsave("Exp_Heatmap_AVG_Res03.pdf", useDingbats = F, height = 10, width = 5, units = "in")