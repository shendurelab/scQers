library(Seurat)
library(MouseGastrulationData)
library(pryr)
library(tidyverse)



# # # # # # # # # # # # # # # # # # # # 
# load the data scQer mEB experiment
# # # # # # # # # # # # # # # # # # # # 

file_mEB_scQer <- "/net/shendure/vol10/projects/JBL/integration_scRNA_mEB_202208/nobackup/gex_scREA_mEB_bottlenecked_all_reps_A_B_2B_subclus_doub_rem_20220819.RDS"
gex_mEB_scQer <- readRDS(file_mEB_scQer)

# add the coarser cluster ids
coarse_clusters <- list(c(0,5),c(1,7),c(4,6),c(2,9,12),3)
names(coarse_clusters) <- c("pluri","ecto","endo","meso","unk")

cluster_ids <- gex_mEB_scQer$seurat_clusters

coarse_cluster_ids <- vector(mode="character",length=length(cluster_ids))
for (coarse_cluster_oi in names(coarse_clusters)){
  cluster_idx <- which(cluster_ids %in% coarse_clusters[[coarse_cluster_oi]])
  coarse_cluster_ids[cluster_idx] <- coarse_cluster_oi
}

gex_mEB_scQer$coarse_cluster_ids <- coarse_cluster_ids

gex_mEB_scQer_minimal <- CreateSeuratObject(counts=gex_mEB_scQer@assays$RNA@counts,
                  assay="RNA", meta.data=gex_mEB_scQer@meta.data %>% 
                                                select(orig.ident,
                                                       nCount_RNA,
                                                       nFeature_RNA,
                                                       percent.mt=percent_mt,
                                                       doublet_score,
                                                       original_coarse_cluster_id=coarse_cluster_ids))

gex_mEB_scQer_minimal$sample <- "mEB_d21_scQer"
rm(gex_mEB_scQer)

saveRDS(gex_mEB_scQer,"minimal_scQer_mEB_d21_20220824.RDS")


print("done loading scQer data")

# # # # # # # # # # # # # # # # # # # # 
# load the data from Pijuan-Sala
# # # # # # # # # # # # # # # # # # # # 

sample_PS_name <- "Pijuan_Sala_2019"
gex_PS <- EmbryoAtlasData(samples = c(1:10,12:20,23:37))

# change gene names
matrix_PS <- gex_PS@assays@data$counts
cells_PS <- Cells(gex_PS)
bool_valid_cells_PS <- !is.na(gex_PS$celltype)

PS_metadata <- data.frame(gex_PS@colData[bool_valid_cells_PS,])
cells_PS2 <- cells_PS[bool_valid_cells_PS]
matrix_PS2 <- matrix_PS[,bool_valid_cells_PS]

# get corresponding gene names
gene_names_mat <- rownames(gex_PS)
feature_table <- read.table("/net/shendure/vol10/projects/JBL/integration_scRNA_mEB_202208/nobackup/features.tsv.gz",header=FALSE,sep="\t")
colnames(feature_table) <- c("ID","gene","type")
feature_table2 <- feature_table %>% select(ID,gene)
rownames(feature_table2) <- feature_table2$ID
bool_all_dups <- function(x){
  return(!(duplicated(x) | duplicated(x, fromLast = TRUE)))
}

gene_names_mat_simple <- feature_table2[gene_names_mat,"gene"]
gene_name_dup <- bool_all_dups(gene_names_mat_simple)
gene_names_mat_simple2 <- gene_names_mat_simple[gene_name_dup]

matrix_PS3 <- matrix_PS2[gene_name_dup,]
rownames(matrix_PS3) <-gene_names_mat_simple2 
colnames(matrix_PS3) <- cells_PS2

# Create a Seurat object with the new gene name
gex_PS_seurat <- CreateSeuratObject(counts=matrix_PS3,meta.data = PS_metadata)
gex_PS_seurat$sample <- sample_PS_name


saveRDS(gex_PS_seurat,"minimal_Pijuan_Sala_20220824.RDS")

print("done loading Pijuan-Sala")


print("Combine objects")

genes_scQer <- rownames(gex_mEB_scQer_minimal)
genes_PS <- rownames(gex_PS_seurat)
genes_common_all <- Reduce(intersect,
                           list(genes_scQer,genes_PS))

gex_scQer2 <- CreateSeuratObject(counts=gex_mEB_scQer_minimal@assays$RNA@counts[genes_common_all,],
                                 meta.data=gex_mEB_scQer_minimal@meta.data)

gex_PS2 <- CreateSeuratObject(counts=gex_PS_seurat@assays$RNA@counts[genes_common_all,],
                                 meta.data=gex_PS_seurat@meta.data)


gex_combined <- merge(gex_scQer2,y=gex_PS2)






