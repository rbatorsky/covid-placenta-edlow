# Analyze DEG, Fig 2 and Fig 3 C,D

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(ggrepel)
  library(openxlsx)
  library(ggpubr)
  
})

source('scripts/00_get_started.R')

# load data ----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE_norm.rds"
so = readRDS(so_file)

# metadata
# Column "sub_cluster_names" has the cell-types
# Column integrated_snn_res.0.3 has the clusters that were used for differential expression
# differential expression was done between the two groups in column "covid_status"

# convert the cluster names
names = so@meta.data %>% 
  dplyr::select(c('integrated_snn_res.0.3','sub_cluster_names')) %>% 
  mutate(cluster = as.character(integrated_snn_res.0.3)) 

names = unique(names)

write.csv(names, "analysis/seurat_object/cluster_names.csv", row.names=F)

# read in DEG  ----

# deg with bhp and no cell number downsampling 
deg = read.csv("analysis/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_fixbh.csv")

deg = deg %>% mutate(cluster = as.character(cluster)) %>% 
  full_join(names, by=c("cluster")) 

deg_filter <- deg %>%
  filter(p_val_adj_bh < 0.05 & abs(avg_log2FC) > 0.2) %>%
  group_by(sub_cluster_names) %>%
  arrange(factor(sub_cluster_names, levels=c("HBC 0","HBC 1","HBC 2","HBC 3","HBC 4","HBC 5","HBC 6","HBC 7","PAMM","Monocyte"))) %>%
  mutate(dir = ifelse(avg_log2FC>0,'up','dn'))

table(deg_filter$sub_cluster_names, deg_filter$dir)

write.xlsx(deg_filter, "analysis/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_fixbh_bh_0.05_lfc_0.2.xlsx")

# how does it change if the bonferroni is used instead, not used
# deg_bon_filter <- deg %>%
#   filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>%
#   group_by(sub_cluster_names) %>%
#   arrange(factor(sub_cluster_names, levels=c("HBC 0","HBC 1","HBC 2","HBC 3","HBC 4","HBC 5","HBC 6","HBC 7","PAMM","Monocyte"))) %>%
#   mutate(dir = ifelse(avg_log2FC>0,'up','dn'))
# 
# table(deg_bon_filter$sub_cluster_names, deg_bon_filter$dir)

um = DimPlot(so, split.by="sample_type", group.by="sub_cluster_names")
ggsave(filename = "analysis/deg/plots/covid_control_umap.pdf",um,width = 8, height = 4)


# DEG and number of cells barplot, Fig 2A ----
# number of deg
tt= data.frame(table(deg_filter$sub_cluster_names, deg_filter$dir))  %>%
  arrange(factor(Var1, levels=c("HBC 0","HBC 1","HBC 2","HBC 3","HBC 4","HBC 5","HBC 6","HBC 7","Monocyte","PAMM")))

tt$Var1 = factor(tt$Var1, levels=c("HBC 0","HBC 1","HBC 2","HBC 3","HBC 4","HBC 5","HBC 6","HBC 7","Monocyte","PAMM"))

tt$cluster = tt$Var1
tt$dir = tt$Var2
tt$Freq = as.numeric(tt$Freq)
tt$Freq = ifelse(tt$dir == "dn",-tt$Freq,tt$Freq)

p = ggplot(data=tt, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="blue", aes(alpha=dir)) +
  scale_alpha_manual(values=c(0.4,1),
                     name="dir",
                     breaks=c("dn","up"),         # can use to set order
                     labels=c("dn","up")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Number of DEG") +
  xlab("")+ theme_classic()+ theme(axis.text=element_text(size=12),
                                                  axis.title=element_text(size=14))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(drop=FALSE) 


ggsave(p, filename = "analysis/deg/plots/ndeg_newnames_updn.pdf",
       device = cairo_pdf,
       width = 4, height = 4,
       units = "in")


# number of cells

p1 = ggplot(data = so@meta.data, aes(x=sub_cluster_names, fill = factor(sample_type))) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values=c('lightgrey','red')) + 
  theme_classic() +
  xlab("") + 
  ylab("Number of DEG") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(legend.position="none", panel.spacing.x = unit(0, "npc"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

print(p1)
ggsave(p1, filename = "analysis/deg/plots/covid_control_fraction.pdf",
       device = cairo_pdf,
       width = 3.3, height = 3.3,
       units = "in")

# GO for COVID recovered vs. COVID negative DEG for Human Placenta, select term dotplot Fig 2B -----

# select term dotplot 
ck = readRDS("analysis/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2.rds")
select_terms_workbook = read.xlsx("analysis/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_with_categories.xlsx")

select_terms_workbook = select_terms_workbook %>%
  dplyr::filter(old_select == 1)

unique(select_terms_workbook$category)
deg_term_select_combine = unique(select_terms_workbook$Description)

ck_select=ck %>% dplyr::filter(Description %in% deg_term_select_combine) %>% 
  dplyr::filter(cluster !="Monocytes")

ck_select@compareClusterResult$Description = factor(ck_select@compareClusterResult$Description, levels=rev(deg_term_select_combine))

select_terms_workbook = select_terms_workbook %>%
  dplyr::filter(old_select == 1) %>%
  mutate(category = factor(category, levels = c("immune", "migration", "metabolism", "stress", "vascular", "ribosome", "protein", "phago", "actin","lipid"))) %>%
  arrange(category)

deg_term_select_combine = unique(select_terms_workbook$Description)

p = ordered_cp_dotplot(ck, deg_term_select_combine)

print(p)
ggsave(p, filename = "analysis/markers/plots/select_terms_go.pdf",
       device = cairo_pdf,
       width = 9, height = 12, 
       units = "in")


# GO bar plot, using combine plots, not used ----
ck_top  = data.frame(ck)  %>% 
  arrange(p.adjust) %>% group_by(cluster) %>% 
  top_n(10, wt=-p.adjust)

ck_top$nlogp = -log10(ck_top$p.adjust)
unique(ck_top$cluster)

color=sample(brewer.pal(8, "Spectral"))
categories=unique(ck_top$Cluster)
names(color) = categories
outdir="analysis/deg/plots/"

p0 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_0"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28) + xlab("")

p1 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_1"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28)+ ylab("")+ xlab("")


p2 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_2"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28)+ ylab("")+ xlab("")

p3 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_3"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28) + ylab("")+ xlab("")

p5 = ggplot(ck_top %>% dplyr::filter(cluster=="HBCs_5"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28) + ylab("")+ xlab("")

ppam = ggplot(ck_top %>% dplyr::filter(cluster=="PAMs"),
              aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic() + 
  ylim(0,28)+ ylab("")+ xlab("")

pmono = ggplot(ck_top %>% dplyr::filter(cluster=="Monocytes"),
            aes(fct_reorder(Description, -nlogp), nlogp, fill = cluster)) + 
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 12)) + 
  coord_flip() + theme_classic()+ 
  ylim(0,28)+ ylab("")+ xlab("")

ggarrange(pmono +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 80) ),
          ppam +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l=110) ), 
          p5 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 20) ), 
          p3 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 110)),
          p2 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l =  40) ),
          p1 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 0)),
          p0 +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(l = 50) ),
          ncol = 1)






# GO BAR plot, another method, not used
ck_top  = data.frame(ck@compareClusterResult)  %>% 
  arrange(p.adjust) %>% group_by(cluster) %>% 
  slice(1:10)

ck_top$nlogp = -log10(ck_top$p.adjust)

color=sample(brewer.pal(8, "Spectral"))
categories=unique(ck_top$Cluster)
names(color) = categories
outdir="analysis/deg/plots/"

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
  scale_y_continuous(expand = c(0,0)) + 
    ylim(0,25)
  
  ggsave(p, filename = paste0(outdir,paste0(i,"_gobar.pdf")),
         device = cairo_pdf,
         width = 9, height = 3,
         units = "in")
  
}


view(ck_top)


# plot individual GO terms for individual heat map, Fig 3 ----

go_terms_select=c("ATP metabolic process",
                  "phagocytosis",
                  "neuroinflammatory response",
                  "neuron death",
                  "regulation of phagocytosis",
                  "cell chemotaxis",
                  "RNA catabolic process",
                  "oxidative phosphorylation",
                  "cellular response to oxidative stress")



term="phagocytosis"
ck = readRDS("/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2.rds")

ck_select = ck %>% dplyr::filter(Description == term)
ck_select_genes = data.frame(genes = ck_select@compareClusterResult$geneID)
ck_select_genes$split = gsub("/",",",ck_select_genes$genes)
ck_select_genes_list = unique(unlist(strsplit(ck_select_genes$split,",")))

## convert the cluster names
names = so@meta.data %>% dplyr::select(c('integrated_snn_res.0.3','sub_cluster_names')) %>% unique() %>% mutate(cluster = as.character(integrated_snn_res.0.3))
rownames(names) = NULL

deg = deg %>% mutate(cluster = as.character(cluster)) %>% 
  full_join(names, by=c("cluster"))  %>%
  mutate(cluster = sub_cluster_names) 

deg_select = deg %>% 
  dplyr::filter(gene %in% ck_select_genes_list) %>%
  arrange(factor(cluster, levels=c("HBCs_0","HBCs_1","HBCs_2","HBCs_3","HBCs_4","HBCs_5","HBCs_6","HBCs_7","PAMs","Monocytes"))) %>%
  dplyr::select("gene","avg_log2FC","p_val_adj_bh", "cluster")

deg_select$cluster = factor(deg_select$cluster, levels=c("HBCs_0","HBCs_1","HBCs_2","HBCs_3","HBCs_4","HBCs_5","HBCs_6","HBCs_7","PAMs","Monocytes"))

lfc_select = deg_select %>% 
  dplyr::select("gene","avg_log2FC", "cluster") %>% 
  spread("cluster","avg_log2FC")

rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL

lfc_select_mat = as.matrix(lfc_select)
lfc_select_mat[is.na(lfc_select_mat)] <- 0

p_select = deg_select %>% 
  dplyr::select("gene","p_val_adj_bh", "cluster") %>% 
  spread("cluster","p_val_adj_bh")
rownames(p_select) = p_select$gene
p_select$gene = NULL
p_select_mat = as.matrix(p_select)
p_select_mat[is.na(p_select_mat)] <- 1



ph = Heatmap(lfc_select_mat, 
             cluster_columns=FALSE, 
             cell_fun = function(j, i, x, y, w, h, fill) {
  if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
    grid.text("*", x, y)
  }
})


show(ph)


pdf("~/human_scrna_edlow_2021/deg/plots/phago_cheat.pdf",
    width = 4.25, height = 10.5)
print(ph)
dev.off()

# DEG GENE HEATMAPS, to show reproduciblity across patients -----
#https://divingintogeneticsandgenomics.rbind.io/post/enhancement-of-scrnaseq-heatmap-using-complexheatmap/

outdir="analysis/deg/plots/"

fc_thresh=1

# avg expression heatmap by donor AND covid status
for (cl in unique(so@meta.data$sub_cluster_names)){
  print(cl)
  cl1 = subset(so, sub_cluster_names == cl)
  cl1@meta.data$donor_status = paste0(cl1@meta.data$covid_status, "_", cl1@meta.data$donor)
  cl1@meta.data$donor_status = factor(cl1@meta.data$donor_status, levels=sort(unique(cl1@meta.data$donor_status)))
  cluster_anno<- cl1@meta.data$donor_status
  cl1.avg <- AverageExpression(cl1, group.by = "donor_status")
  counts = data.frame(cl1.avg$RNA)
  counts$gene = rownames(counts)
  
  cl1_top_genes <- deg_filter %>% 
    dplyr::filter(sub_cluster_names == cl) %>% 
    dplyr::filter(abs(avg_log2FC) > fc_thresh) %>%
    top_n(n = 25, wt = -p_val_adj) %>% dplyr::filter(gene != "XIST")
  
  counts = counts %>% dplyr::filter(gene %in% cl1_top_genes$gene)
  counts$gene=NULL
  p = pheatmap(counts, scale="row", cluster_cols = T)
  
  pdf(paste0(outdir,cl,"_top25_absg",fc_thresh,"_donorstat_pheatmap.pdf"),width = 7, height = 5)
  print(p)
  dev.off()


}

# single cell heatmap by Donor and COVID status only 
for (cl in unique(so@meta.data$sub_cluster_names)){
  print(cl)
  cl1 = subset(so, sub_cluster_names == cl)
  cl1@meta.data$donor_status = paste0(cl1@meta.data$covid_status, "_", cl1@meta.data$donor)
  cl1@meta.data$donor_status = factor(cl1@meta.data$donor_status, levels=sort(unique(cl1@meta.data$donor_status)))
  cluster_anno<- cl1@meta.data$donor_status
  
  cl1_top_genes <- deg_filter %>% 
    dplyr::filter(sub_cluster_names == cl) %>% 
    dplyr::filter(abs(avg_log2FC) > fc_thresh) %>%
    top_n(n = 25, wt = -p_val_adj) %>% dplyr::filter(gene != "XIST")
  
  mat<- cl1[["RNA"]]@data[cl1_top_genes$gene, ] %>% as.matrix()
  
  ## scale the rows
  mat<- t(scale(t(mat)))
  
  ## annotation
  quantile(mat, c(0.1, 0.95))
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c"))
  
  p = Heatmap(mat, name = "Expression",  
              column_split = factor(cluster_anno),
              cluster_columns = F,
              show_column_dend = FALSE,
              cluster_column_slices = F,
              column_title_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = T,
              show_row_dend = T,
              col = col_fun,
              row_names_gp = gpar(fontsize = 8),
              column_title_rot = 90,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(12)))),
              show_column_names = FALSE,
              use_raster = F) 
  
  
  
  pdf(paste0(outdir,cl,"_top25_absg",fc_thresh,"_donorstat_heatmap.pdf"),width = 7, height = 5)
  print(p)
  dev.off()
  
  
}


# single cell heatmap by COVID status only
for (cl in unique(so@meta.data$sub_cluster_names)){
  print(cl)
  cl1 = subset(so, sub_cluster_names == cl)
  cl1.avg <- AverageExpression(cl1, group.by = "covid_status")
  counts = data.frame(cl1.avg$RNA)
  counts$gene = rownames(counts)
  
  cl1_top_genes <- deg_filter %>% 
    dplyr::filter(sub_cluster_names == cl) %>% 
    dplyr::filter(abs(avg_log2FC) > fc_thresh) %>%
    top_n(n = 25, wt = -p_val_adj) %>% dplyr::filter(gene != "XIST")
  
  mat<- cl1[["RNA"]]@data[cl1_top_genes$gene, ] %>% as.matrix()
  
  ## scale the rows
  mat<- t(scale(t(mat)))
  
  ## annotation
  cluster_anno<- cl1@meta.data$covid_status
  
  #col_fun = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c") 
  #circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  quantile(mat, c(0.1, 0.95))
  #col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c"))
  
  p = Heatmap(mat, name = "Expression",  
              column_split = factor(cluster_anno),
              cluster_columns = F,
              show_column_dend = FALSE,
              cluster_column_slices = F,
              column_title_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = T,
              show_row_dend = T,
              col = col_fun,
              row_names_gp = gpar(fontsize = 8),
              column_title_rot = 0,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("grey","red")))),
              show_column_names = FALSE,
              use_raster = F) 

  pdf(paste0(outdir,cl,"_top25_absg",fc_thresh,"_covidstat_heatmap.pdf"),width = 7, height = 5)
  print(p)
  dev.off()
  
  
}


# QC of  DEG results ---- 
# 
# ## No mincell
# no_mincell_fixbh_results = read.csv("/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_fixbh.csv")
# head(no_mincell_fixbh_results)
# 
# # rawp
# sig_results_rawp = no_mincell_fixbh_results %>% dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_rawp$cluster)
# 
# # bhp
# sig_results_fix_bh = no_mincell_fixbh_results %>% dplyr::filter(p_val_adj_bh < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_fix_bh$cluster)
# #write.csv(sig_results_fix_bh, "~/slonimlab/rbator01/human_scrna_edlow_2021/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_bh_0.05_lfc_0.2_fixbh.csv")
# 
# # bonp
# sig_results_bon = no_mincell_fixbh_results %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_bon$cluster)
# 
# 
# ## Yes Mincell
# fix_bh_results = read.csv("~/slonimlab/rbator01/human_scrna_edlow_2021/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_mincell_fixbh.csv")
# 
# #rawp
# sig_results_rawp = fix_bh_results %>% dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_rawp$cluster)
# 
# # bhp
# sig_results_fix_bh = fix_bh_results %>% dplyr::filter(p_val_adj_bh < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_fix_bh$cluster)
# 
# #bonp
# sig_results_bon = fix_bh_results %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% arrange(cluster)
# table(sig_results_bon$cluster)
# 

## Checking some of Britt's markers
# wilcox_markers=read.csv("~/slonimlab/rbator01/human_scrna_edlow_2021/markers.wilcox_seuratobj_subset_integrated_HBCs.csv",row.names=NULL)
# negbin_markers=read.csv("~/slonimlab/rbator01/human_scrna_edlow_2021/markers.negbin_seuratobj_subset_integrated_HBCs.csv")
# wilcox_markers_prenorm=read.csv("~/slonimlab/rbator01/human_scrna_edlow_2021/before_norm/markers.wilcox_seuratobj_subset_integrated_HBCs.csv")
# negbin_markers_prenorm=read.csv("~/slonimlab/rbator01/human_scrna_edlow_2021/before_norm/markers.negbin_seuratobj_subset_integrated_HBCs.csv")
# britt_wilcox=read.csv("/cluster/home/rbator01/slonimlab/rbator01/human_scrna_edlow_2021/markers_res03_integrated_wilcox.csv")
# britt_negbin=read.csv("/cluster/home/rbator01/slonimlab/rbator01/human_scrna_edlow_2021/markers_res03_integrated_negbinom.csv")


# adding in 7/15: There are two DEG files, we want to perform go analysis and heatmaps on both of them and decide which is the better set

# This on uses downsampling to the minimum cells in either group (covid rec, covid neg) - remove this
# deg_ds = read.csv("/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_mincell_fixbh.csv")
# 
# deg_ds = deg_ds %>% mutate(cluster = as.character(cluster)) %>%
#   full_join(names, by=c("cluster"))
# 
# head(deg_ds)
# deg_ds_filter <- deg_ds %>%
#   filter(p_val_adj_bh < 0.05 & abs(avg_log2FC) > 0.2) %>%
#   group_by(sub_cluster_names) %>%
#   arrange(factor(sub_cluster_names, levels=c("HBC 0","HBC 1","HBC 2","HBC 3","HBC 4","HBC 5","HBC 6","HBC 7","PAMM","Monocyte")))
# 
# table(deg_ds_filter$sub_cluster_names)


# mark manually selected terms with category
# select_terms_workbook = read.xlsx("analysis/deg/AGE categories ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2.xlsx")
# 
# #deg_term_select=c(imm, migration, meta, stress, vascular, ribo, prot, phago, actin,lipid)
# 
# ck_df = ck_df %>%
#   mutate(old_select = ifelse(Description %in% deg_term_select_combine, 1,0)) %>%
#   mutate(category = ifelse(Description %in% strsplit(imm,'\n')[[1]], 'immune',
#                            ifelse(Description %in% strsplit(migration,'\n')[[1]], 'migration',
#                                   ifelse(Description %in% strsplit(meta,'\n')[[1]], 'metabolism',
#                                          ifelse(Description %in% strsplit(stress,'\n')[[1]], 'stress',
#                                                 ifelse(Description %in% strsplit(vascular,'\n')[[1]],'vascular',
#                                                        ifelse(Description %in% strsplit(ribo,'\n')[[1]],'ribosome',
#                                                               ifelse(Description %in% strsplit(prot,'\n')[[1]],'protein',
#                                                                      ifelse(Description %in% strsplit(phago,'\n')[[1]],'phago',
#                                                                             ifelse(Description %in% strsplit(actin,'\n')[[1]],'actin',
#                                                                                    ifelse(Description %in% strsplit(lipid,'\n')[[1]], 'lipid','NA')))))))))))
# 
# 
# 
# 
# 
# view(ck_df %>% dplyr::filter(old_select == 1 ))
# write.xlsx(ck_df, "analysis/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_with_categories.xlsx")

