# This script plots expression of several gene modules from published data, Fig 3E

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

source('scripts/00_get_started.R')

# read in data ----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE.rds"

gene.sets <- read.csv("module/02_Signatures_HumanConverted_Collated.csv")

# Askenase, Pulido, Kracht paper modules ----
Askenase_CD14monocytes <- list(gene.sets$Askenase_CD14monocytes)
Askenase_microglia <- list(gene.sets$Askenase_microglia)
Askenase_BorderAssociatedMacrophages <- list(gene.sets$Askenase_BorderAssociatedMacrophages)
Pulido_up_treated_microglia_Fig1 <- list(gene.sets$Pulido_up_treated_microglia_Fig1)
Pulido_pos_FC_up_LPS <- list(gene.sets$Pulido_pos_FC_up_LPS)
Kracht_microglia_activation_fig2C <- list(gene.sets$Kracht_microglia_activation_fig2C)

so <- AddModuleScore(object = so, features = Askenase_microglia, name = "Askenase_microglia", seed.use = 1988)
so <- AddModuleScore(object = so, features = Askenase_CD14monocytes, name = "Askenase_CD14monocytes", seed.use = 1988)
so <- AddModuleScore(object = so, features = Askenase_BorderAssociatedMacrophages, name = "Askenase_BorderAssociatedMacrophages", seed.use = 1988)
so <- AddModuleScore(object = so, 
                     features = Pulido_up_treated_microglia_Fig1, 
                     name = "Pulido_up_treated_microglia_Fig1", seed.use = 1988)
so <- AddModuleScore(object = so, 
                     features = Pulido_pos_FC_up_LPS, 
                     name = "Pulido_pos_FC_up_LPS", seed.use = 1988)
so <- AddModuleScore(object = so, 
                     features = Kracht_microglia_activation_fig2C, name = "Kracht_microglia_activation_fig2C", seed.use = 1988)
colnames(so@meta.data)
Idents(so) = "sub_cluster_names"


# Thomas paper module scores ---
markers.bao <- read.xlsx("~/slonimlab/rbator01/reference_data/thomas_2020/jem_20200891_Sdat_combined.xlsx")

clusters=unique(markers.bao$cluster)
view(markers.bao)
for (cl in clusters){
  print(cl)
  cluster_filt = markers.bao %>% dplyr::filter(cluster == cl)
  feature=cluster_filt$Gene
  
  all_genes=rownames(so@assays$RNA@counts)
  
  feature=intersect(feature,all_genes)
  feature=list(c(feature))
  so <- AddModuleScore(
    object = so,
    features = feature,
    name = cl
  )
}

features = c("YS1","Mono1","CS101","SAMac1","Kup_cell1" )

cols.HBCs <- c("#60dca0",
               "#4f91bf",
               "#81d94c",
               "#76d5d3",
               "#3b5b2c",
               "#a9cd84",
               "#d65232",
               "#b6475a",      
               "#519279",
               "#51943e")

DefaultAssay(so) = "RNA"


# don't use violin plots, confusing
# VlnPlot(so, features = c("Askenase_microglia1"),
#         pt.size = 0, cols = cols.HBCs,
#         slot="data",
#         assay="RNA")
# 
# 
# VlnPlot(so, features = c("Askenase_CD14monocytes1"),
#         pt.size = 0, cols = cols.HBCs,
#         slot="data",
#         assay="RNA")
# 
# 
# VlnPlot(so, features = c("Askenase_BorderAssociatedMacrophages1"),
#         pt.size = 0, cols = cols.HBCs,
#         slot="data",
#         assay="RNA")


# show as a heatmap

features = c("Askenase_microglia1","Askenase_BorderAssociatedMacrophages1", "YS1","Mono1")

so[['module']]= CreateAssayObject(data = t(x = FetchData(object = so, vars = features)))
so_avg <- AverageExpression(so)
full=t(so_avg$module)
head(full)

scale_mod = t(scale(t(so_avg$module)))

ph = Heatmap(scale_mod, 
             cluster_columns=T,
             color = rev(brewer.pal(n = 7, name = "RdYlBu"))
)



pdf("module/4mod_cheat.pdf",
    width = 7, height = 3)
print(ph)
dev.off()


p=pheatmap(so_avg$module,
           scale='row',
           cluster_cols = F,
           cluster_rows = T) 

pdf("module/pheatmap_3mod.pdf", 
    height=2.5, width = 5)
print(p)
dev.off()


# Microglial activation modules -----
Pulido_up_treated_microglia_Fig1 <- list(gene.sets$Pulido_up_treated_microglia_Fig1)
so <- AddModuleScore(object = so, 
                     features = Pulido_up_treated_microglia_Fig1, 
                     name = "Pulido_up_treated_microglia_Fig1", seed.use = 1988)

VlnPlot(so, features = c("Pulido_up_treated_microglia_Fig11"),
        pt.size = 0, cols = cols.HBCs)


so@meta.data$donor_status = paste0(so@meta.data$covid_status, "_", so@meta.data$donor)
so@meta.data$donor_status = factor(so@meta.data$donor_status, levels=sort(unique(so@meta.data$donor_status)))

colnames(so@meta.data)
Idents(so) = "sub_cluster_names"

VlnPlot(so, features = c("Pulido_up_treated_microglia_Fig11"),
        pt.size = 0)

cols.status <- c("grey","red")

VlnPlot(so, features = c("Pulido_up_treated_microglia_Fig11"),
        pt.size = 0, split.by = "covid_status", cols = cols.status)


