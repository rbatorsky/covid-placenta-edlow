# This script reads in the seurat object, checks normalization and makes summary HC and PCA plots

suppressPackageStartupMessages({library(org.Hs.eg.db)
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(Libra)
  library(DESeq2)
  library(magrittr)
  library(ggrepel)
  library("pheatmap")
  library("RColorBrewer")
  })

source('scripts/00_get_started.R')

# Load integrated dataset from Britt ----
so_file="analysis/seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE.rds"
so = readRDS(so_file)
table(so@meta.data$donor, so@meta.data$covid_status)

res="sub_cluster_names"
Idents(so) = res

# Since these are equal, this has not been normalized
all(so[['RNA']]@counts@x == so[['RNA']]@data@x)

# Normalize
DefaultAssay(so) <- "RNA"
so = NormalizeData(so)

# Rename the clusters 
# convert the cluster names
names = so@meta.data %>% 
  dplyr::select(c('integrated_snn_res.0.3','sub_cluster_names')) %>% 
  mutate(cluster = as.character(integrated_snn_res.0.3)) 

rownames(names) = rownames(so@meta.data)

names$sub_cluster_names = gsub(pattern = "HBCs_", replacement = "HBC ",x = names$sub_cluster_names)
names$sub_cluster_names = gsub(pattern = "Monocytes", replacement = "Monocyte",x = names$sub_cluster_names)
names$sub_cluster_names = gsub(pattern = "PAMs", replacement = "PAMM",x = names$sub_cluster_names)

so = AddMetaData(so,
                 names$sub_cluster_names,
                 "sub_cluster_names")

saveRDS(so, "seurat_object/seuratobj_subset_integrated_HBCs_USETHISONE_norm.rds")

# PCA for quality control ----

# assign one cluster to do full sample pseudobulk
so@meta.data$one_cluster = "one"

unique(so@meta.data$donor)
unique(so@meta.data$covid_status)

level1='negative'
level2='recovered'

metadata <- so@meta.data
metadata$donor= as.factor(metadata$donor)
metadata$covid_status = as.factor(metadata$covid_status)
metadata$sub_cluster_names = as.factor(metadata$sub_cluster_names)
so@meta.data = metadata

min_cells = 0
min_reps = 0
cell_type_col="one_cluster"
replicate_col="donor"
label_col="covid_status"
level1=level1
level2=level2

pseudobulks_orig_ident=Libra::to_pseudobulk(so, 
                                            min_cells = min_cells,
                                            min_reps = min_reps, 
                                            cell_type_col=cell_type_col,
                                            replicate_col=replicate_col,
                                            label_col=label_col)

rownames(metadata) = NULL
deseq_meta = unique(so@meta.data[c("donor", "sex","covid_status","batch")])
deseq_meta
rownames(deseq_meta) = deseq_meta$donor

deseq_counts = pseudobulks_orig_ident$one
colnames(deseq_counts)

colnames = data.frame(names=colnames(deseq_counts))
colnames = colnames %>% separate(names, into=c("name","status"), sep=":")
colnames(deseq_counts) = colnames$name
deseq_counts <- deseq_counts[, rownames(deseq_meta)]
colnames(deseq_counts)

# check order
all(colnames(deseq_counts) %in% rownames(deseq_meta))
all(colnames(deseq_counts) == rownames(deseq_meta))

## run deseq analysis
#design = model.matrix(~ group, data = targets)

dds = DESeqDataSetFromMatrix(countData = deseq_counts,
                             colData = deseq_meta,
                             design = ~ sex)


dds = DESeq(dds)

# Plotting a gene
d <- plotCounts(dds, gene="CCL8", intgroup="sex", returnData=TRUE)

ggplot(d, aes(x = sex, y = count, color = sex)) +
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) +
  theme_bw() +
  ggtitle("Hnrnpm") +
  theme(plot.title = element_text(hjust = 0.5))


res = results(dds)
res.df = data.frame(res)
res.df$gene = rownames(res.df)
res.df %>% dplyr::filter(gene == "CCL8")

# check the PCA
vsd <- DESeq2::vst(dds)
rv=rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(100, length(rv)))]
pc = prcomp(t(assay(vsd)[select,]))
loadings = as.data.frame(pc$rotation)
aload = abs(loadings)
sweep(aload, 2, colSums(aload), "/")
aload$gene = rownames(aload)
 
aload_long <- aload %>% 
  gather(pc, value, -gene) %>%
  group_by(pc) %>%
  top_n(20, wt=value)

df <- as.data.frame(colData(dds)[,c("donor","sex","covid_status", "batch")])


annoCol<-list(sex=c(M="blue", F="red"),
              covid_status=c(negative="orange", recovered="green"),
              batch=c(AEb14="blue", AEb13 = "yellow", AEb12="grey", AEb11="green",AEb9 = "orange",
                      AEb8 = "black", AEb7 = "purple"))


pheatmap(assay(vsd)[select,], 
         cluster_rows=T, 
         show_rownames=T,
         cluster_cols=T, 
         annotation_col=df, 
         annotation_colors = annoCol)

vsd_cor <- cor(assay(vsd))
pheatmap(vsd_cor, 
         cluster_rows=T, 
         show_rownames=T,
         cluster_cols=T, 
         annotation_col=df, 
         annotation_colors = annoCol)

# plotPCA ----
vsd <- DESeq2::vst(dds)
rv=rowVars(assay(vsd))
intgroup=c("donor","covid_status")
pcaData <- plotPCA(vsd, intgroup=intgroup, returnData=TRUE)
ntop=10
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])

group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = " : "))
}else{
  colData(object)[[intgroup]]
}

d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, name = colData(vsd)[,1])


attr(d, "percentVar") <- percentVar[1:2]


pcaData <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                      group = group, intgroup.df, name = colnames(vsd))


head(pcaData)
percentVar <- round(100 * percentVar)

ggplot(pcaData, aes(PC1, PC2, color=covid_status,label=donor)) +
  geom_point(size=3) +
  geom_text(hjust=0.2, vjust=0.2) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


ggplot(pcaData, aes(PC3, PC4, color=covid_status,label=donor)) +
  geom_point(size=3) +
  geom_text(hjust=0.2, vjust=0.2) + 
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
  coord_fixed()
