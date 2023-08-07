suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(openxlsx)
})

source('scripts/00_get_started.R')

# load data -----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE.rds"
so = readRDS(so_file)
res="sub_cluster_names"
Idents(so) = res
DefaultAssay(so) <- "RNA"

clusters=unique(so@meta.data$sub_cluster_names)
colnames(so@meta.data)

# calc markers for each cluster ----
results=NULL

for (cl in clusters){
  de_cl <- FindMarkers(so, ident.1 = cl, only.pos = TRUE, logfc.threshold = 0)
  de_cl$cluster = cl
  de_cl$gene = rownames(de_cl)
  de_cl = de_cl %>% dplyr::select(c('gene',everything()))
  
  
  if (is.null(results)){
    results = de_cl
  }else{
    results = rbind(results, de_cl)
  }
}

write.xlsx(results,"markers/markers_seuratobj_subset_integrated_HBCs_USETHISONE_norm_lfc_0.xlsx")

# hbc vs. pamm markers ----
so_file="seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE_norm.rds"
so = readRDS(so_file)
res="sub_cluster_names"
Idents(so) = res
DefaultAssay(so) <- "RNA"

markers_hbc_pamm <- FindMarkers(object = so,
                                ident.1 = c("HBCs_1","HBCs_2","HBCs_0","HBCs_3","HBCs_4","HBCs_6","HBCs_5","HBCs_7"),
                                ident.2 = c('PAMs'))

p_thresh=0.05
lfc_thresh=1

markers_hbc_pamm$gene = rownames(markers_hbc_pamm)

hbc_pamm_sig = markers_hbc_pamm %>%
  dplyr::filter(p_val_adj<p_thresh & abs(avg_log2FC) > lfc_thresh) %>%
  mutate(direction = ifelse(hbc_pamm_sig$avg_log2FC >0, "HBC","PAMM")) %>%
  dplyr::select(c(gene,direction, everything()))

write.xlsx(hbc_pamm_sig, "markers/seuratobj_subset_integrated_HBCs_USETHISONE_norm_hbc_vs_pamm_markers.xlsx")

# #calc markers WITH subsampling min number of cells per cluster, not used ----
# 
# results=NULL
# for (cl in clusters){
#   so.subset = subset(so, idents=cl)
#   Idents(so.subset)<-"covid_status"
#   min_cell = min(table(so.subset$covid_status))
#   de_cl <- FindMarkers(so.subset, ident.1 = "recovered", ident.2 = "negative",
#                        test.use="MAST", latent.vars = c("donor"), logfc.threshold=0,
#                        max.cells.per.ident = min_cell)
#   de_cl$cluster = cl
#   de_cl$gene = rownames(de_cl)
#   de_cl = de_cl %>% dplyr::select(c('gene',everything()))
#   print(cl)
#   print(length(de_cl$p_val))
#   
#   if(length(de_cl$p_val)>1){
#     de_cl$p_val_adj_bh = p.adjust(de_cl$p_val, method = "BH")
#   }else{
#     de_cl$p_val_adj_bh = de_cl$p_val
#   }
#   
#   if (is.null(results)){
#     results = de_cl
#   }else{
#     results = rbind(results, de_cl)
#   }
# }
# 
# write.xlsx(results,"deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_mincell.xlsx")
# 
# 
# # calc markers WITHOUT subsampling min number of cells per cluster, not used ----
# 
# results=NULL
# 
# for (cl in clusters){
#   so.subset = subset(so, idents=cl)
#   Idents(so.subset)<-"covid_status"
#   min_cell = min(table(so.subset$covid_status))
#   de_cl <- FindMarkers(so.subset, ident.1 = "recovered", ident.2 = "negative",
#                        test.use="MAST", latent.vars = c("donor"), logfc.threshold=0)
#   de_cl$cluster = cl
#   de_cl$gene = rownames(de_cl)
#   de_cl = de_cl %>% dplyr::select(c('gene',everything()))
#   print(cl)
#   print(length(de_cl$p_val))
#   
#   if(length(de_cl$p_val)>1){
#     de_cl$p_val_adj_bh = p.adjust(de_cl$p_val, method = "BH")
#   }else{
#     de_cl$p_val_adj_bh = de_cl$p_val
#   }
#   
#   if (is.null(results)){
#     results = de_cl
#   }else{
#     results = rbind(results, de_cl)
#     
#   }
# }
# 
# write.xlsx(results,"deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE.xlsx")
# 
# 
# # calc markers WITHOUT subsampling min number of cells per cluster - REMOVE FEMALE 171 control -----
# 
# so_no171 = subset(so, subset = donor == "D171", invert =TRUE)
# clusters=unique(so_no171@meta.data$sub_cluster_names)
# unique(so_no171@meta.data$donor)
# results=NULL
# 
# for (cl in clusters){
#   so.subset = subset(so_no171, idents=cl)
#   Idents(so.subset)<-"covid_status"
#   de_cl <- FindMarkers(so.subset, ident.1 = "recovered", ident.2 = "negative",
#                        test.use="MAST", latent.vars = c("donor"), logfc.threshold=0)
#   
#   view(de_cl)
#   de_cl$cluster = cl
#   de_cl$gene = rownames(de_cl)
#   de_cl = de_cl %>% dplyr::select(c('gene',everything()))
#   print(length(de_cl$p_val))
#   
#   if(length(de_cl$p_val)>1){
#     de_cl$p_val_adj_bh = p.adjust(de_cl$p_val, method = "BH")
#   }else{
#     de_cl$p_val_adj_bh = de_cl$p_val
#   }
#   
#   if (is.null(results)){
#     results = de_cl
#   }else{
#     results = rbind(results, de_cl)
#     
#   }
# }
# 
# write.xlsx(results,"deg/deg_rec_vs_neg_lfcthersh0_seuratobj_subset_integrated_HBCs_USETHISONE_no171.xlsx")
# 