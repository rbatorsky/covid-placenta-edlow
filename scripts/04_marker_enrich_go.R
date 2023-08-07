# calculate GO enrichment on cluster markers, cluster DEG, and HBC vs. PAMM markers
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(openxlsx)
})

source('scripts/00_get_started.R')

# marker gene enrich go, with log fold change cutoff ----
markers = read.csv("analysis/markers/markers_seuratobj_subset_integrated_HBCs_USETHISONE.csv")

markers_filter = markers %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.25)
table(markers_filter$cluster)

# convert to ensembl
convert_ens =read.csv("data/AEb12/cellranger/AEb12_CTR/outs/raw_feature_bc_matrix/features.tsv.gz",
                      sep="\t",
                      col.names=c("ens","gene","na"),
                      header=F)

convert_ens = convert_ens %>% dplyr::select(-c('na'))

# It has some duplicated genes
subset(convert_ens,duplicated(gene))
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

# so we make it unique, which is the same process that was done on the original data
convert_ens$gene = make.unique(convert_ens$gene)
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

markers_filter = markers_filter %>% dplyr::left_join(convert_ens,by="gene")

# check that everything converted
markers_filter %>% dplyr::filter(is.na(ens))

table(markers_filter$cluster)

# 3. run clusterprofiler on ensembl ID's

ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = markers_filter,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP",
                     universe=convert_ens$gene
)


saveRDS(ck, "analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25_redo.rds")
write.csv(ck, "analysis/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25_redo.csv")

#ck = readRDS("analysis/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25_redo.rds")
#ck@compareClusterResult %>% dplyr::filter(cluster == "HBCs_2")


p = clusterProfiler::dotplot(ck, x =~cluster,showCategory =10, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

show(p)
ggsave(p, filename = "analysis/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_top20.pdf",
       device = cairo_pdf,
       width = 13, height = 18,
       units = "in")


# marker gene enrich go with NO log fold change cutoff----
markers = read.xlsx("analysis/markers/markers_seuratobj_subset_integrated_HBCs_USETHISONE_norm_lfc_0.xlsx")

markers_filter = markers %>% 
  dplyr::filter(p_val_adj<0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 200)
  

# convert to ensembl
convert_ens =read.csv("data/AEb12/cellranger/AEb12_CTR/outs/raw_feature_bc_matrix/features.tsv.gz",
                      sep="\t",
                      col.names=c("ens","gene","na"),
                      header=F)

convert_ens = convert_ens %>% dplyr::select(-c('na'))

# It has some duplicated genes
subset(convert_ens,duplicated(gene))
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

# so we make it unique, which is the same process that was done on the original data
convert_ens$gene = make.unique(convert_ens$gene)
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

markers_filter = markers_filter %>% dplyr::left_join(convert_ens,by="gene")

# check that everything converted
markers_filter %>% dplyr::filter(is.na(ens))


table(markers_filter$cluster)

# 3. run clusterprofiler on ensembl ID's

ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = markers_filter,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP",
                     universe=convert_ens$gene
)


saveRDS(ck, "analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_top200_lfc_0.rds")
write.csv(ck, "analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_top200_lfc_0.csv")


ck = readRDS("analysis/markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_top200_lfc_0.rds")

# unselected
p = clusterProfiler::dotplot(ck, x =~cluster,showCategory =10, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

show(p)
# ggsave(p, filename = "analysis/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_top20.pdf",
#        device = cairo_pdf,
#        width = 13, height = 18,
#        units = "in")




# deg enrich go ----

so = readRDS("analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE_norm.rds")

# This on uses downsampling to the minimum cells in either group (covid rec, covid neg)
deg = read.csv("analysis/deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_fixbh.csv")

## convert the gene names
names = so@meta.data %>% dplyr::select(c('integrated_snn_res.0.3','sub_cluster_names')) %>%
  unique() %>%
  mutate(cluster = as.character(integrated_snn_res.0.3))

rownames(names) = NULL

deg = deg %>% mutate(cluster = as.character(cluster)) %>%
  full_join(names, by=c("cluster"))  %>%
  mutate(cluster = sub_cluster_names)

## more stringent p_val_adj
# deg_filter <- deg %>%
#   filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>%
#   group_by(cluster)


table(deg_filter$cluster)
deg_filter <- deg %>%
  filter(p_val_adj_bh < 0.05 & abs(avg_log2FC) > 0.2) %>%
  group_by(cluster) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::filter(count>10)

table(deg_filter$cluster)

# No downsample
#0   1   2   3   4   5   6   7   9
#357 723 187 122   3 566  67 177   2


###########
### COVID recovered vs. COVID negative DEG for Human Placenta
##########


# 2. Convert Symbols to ensembl id's using the original cellranger alignment

## Read in features in order to use ensembl IDs
convert_ens =read.csv("/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/raw_data/AEb12/cellranger/AEb12_CTR/outs/raw_feature_bc_matrix/features.tsv.gz",
                      sep="\t",
                      col.names=c("ens","gene","na"),
                      header=F)

convert_ens = convert_ens %>% dplyr::select(-c('na'))

## It has some duplicated genes
subset(convert_ens,duplicated(gene))
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

## so we make it unique, which is the same process that was done on the original data
convert_ens$gene = make.unique(convert_ens$gene)
convert_ens %>% dplyr::filter(grepl('TBCE',gene))

deg_filter = deg_filter %>% dplyr::left_join(convert_ens,by="gene")

# check that everything converted
deg_filter %>% dplyr::filter(is.na(ens))


table(deg_filter$cluster)
# 3. run clusterprofiler on ensembl ID's

ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = deg_filter,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP",
                     universe=convert_ens$gene
)


# # 4. write to file as rds and csv
### Note that the one I gave to andrea was using ENS and no universe.
# ck <- compareCluster(geneCluster = ens ~ cluster,
#                      data = deg_filter,
#                      OrgDb = org.Hs.eg.db,
#                      keyType="ENSEMBL",
#                      fun = "enrichGO",
#                      ont="BP",
#                      readable=TRUE
# )
#


saveRDS(ck, "/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_mincell_fixbh_padj_0.05_lfc_0.2_symbol.rds")
write.csv(ck, "/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_mincell_fixbh_padj_0.05_lfc_0.2_symbol.csv")


p = clusterProfiler::dotplot(ck, x =~cluster,showCategory =5, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

show(p)
ggsave(p, filename = "/cluster/tufts/slonimlab/rbator01/human_scrna_edlow_2021/deg/ck_deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_nomincell_fixbh_padj_0.05_lfc_0.2_top20.pdf",
       device = cairo_pdf,
       width = 13, height = 18,
       units = "in")


# PAMM vs. HBC markers ----

hbc_pamm_sig = read.xlsx("markers/seuratobj_subset_integrated_HBCs_USETHISONE_norm_hbc_vs_pamm_markers.xlsx")


# GO enrichment ----
ck_hbc_pamm<- compareCluster(geneCluster = gene~direction,
                             data = hbc_pamm_sig,
                             OrgDb = org.Hs.eg.db,
                             keyType="SYMBOL",
                             fun = "enrichGO",
                             ont="BP")

saveRDS(ck_hbc_pamm,
        file = "markers/seuratobj_subset_integrated_HBCs_USETHISONE_norm_hbc_vs_pamm_markers_go_bp.rds")
write.csv(ck_hbc_pamm,
          file = "markers/seuratobj_subset_integrated_HBCs_USETHISONE_norm_hbc_vs_pamm_markers_go_bp.xlsx",
          row.names=F)

p = clusterProfiler::dotplot(ck_hbc_pamm,showCategory =10, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# make plots  -----
pamm_markers=c("APOE",
               "APOC1",
               "VIM",
               "G0S2",
               "LGALS1",
               "LYZ",
               "S100A10",
               "S100A6",
               "LGALS3",
               "GPNMB",
               "HLA-DPB1")

hbc_markers=c("LYVE1",
              "F13A1",
              "SELENOP",
              "MRC1",
              "DAB2",
              "CCL2",
              "PLTP",
              "FOLR2",
              "MAF",
              "RNASE1",
              "CCL4")

top_marker_genes=c(pamm_markers, hbc_markers)

DotPlot(so, features = top_marker_genes) +
  xlab("Genes") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


so@meta.data$hbc_pamm = ifelse(so@meta.data$sub_cluster_names %in% c("HBCs_1","HBCs_2","HBCs_0","HBCs_3","HBCs_4","HBCs_6","HBCs_7","HBCs_5"), "HBC",
                               ifelse(so@meta.data$sub_cluster_names %in% c("PAMs"), 'PAMs',
                                      ifelse(so@meta.data$sub_cluster_names %in% c("Monocytes"), 'Monocytes','NA')))

table(hbc_pamm)

Idents(so) = "hbc_pamm"

DotPlot(so, features = top_marker_genes) +
  xlab("Genes") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

thomas_markers = c("FOLR2","HLA-DPB1","CD163","CCR2","CD9","CD68")

DotPlot(so, features = thomas_markers) +
  xlab("Genes") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


FeaturePlot(so, features = thomas_markers)



# label manually selected clusters 
# marker_term_select=c(imm, meta, stress, vascular, ribo, prot, migration, phago, act)
# all=c()
# for (i in marker_term_select){
#   t = strsplit(i, "\n")
#   print(t)
#   all = c(all, t[[1]])
# }
#
# marker_term_select_combine = all
#
# ck_df = ck_df %>%
#   mutate(old_select = ifelse(Description %in% marker_term_select_combine, 1,0)) %>%
#   mutate(category = ifelse(Description %in% strsplit(imm,'\n')[[1]], 'immune',
#                            ifelse(Description %in% strsplit(meta,'\n')[[1]], 'metabolic',
#                                   ifelse(Description %in% strsplit(stress,'\n')[[1]], 'stress',
#                                          ifelse(Description %in% strsplit(ribo,'\n')[[1]], 'ribo',
#                                                 ifelse(Description %in% strsplit(prot,'\n')[[1]],'protein folding',
#                                                        ifelse(Description %in% strsplit(phago,'\n')[[1]],'phago',
#                                                               ifelse(Description %in% strsplit(migration,'\n')[[1]],'migration',
#                                                                      ifelse(Description %in% strsplit(vascular,'\n')[[1]],'vascular',
#                                                                             ifelse(Description %in% strsplit(act,'\n')[[1]], 'actin','NA'))))))))))
#
#
#
#
#
# view(ck_df %>% dplyr::filter(old_select == 1 ))
# write.xlsx(ck_df, "markers/ck_markers_seuratobj_subset_integrated_HBCs_USETHISONE_padj_0.05_lfc_0.25.xlsx")
#
#

