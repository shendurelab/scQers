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


# # # # # # # # # # # # # # # # # # 
# Combine objects
# # # # # # # # # # # # # # # # # # 

genes_scQer <- rownames(gex_mEB_scQer_minimal)
genes_PS <- rownames(gex_PS_seurat)
genes_common_all <- Reduce(intersect,
                           list(genes_scQer,genes_PS))

gex_scQer2 <- CreateSeuratObject(counts=gex_mEB_scQer_minimal@assays$RNA@counts[genes_common_all,],
                                 meta.data=gex_mEB_scQer_minimal@meta.data)

gex_PS2 <- CreateSeuratObject(counts=gex_PS_seurat@assays$RNA@counts[genes_common_all,],
                              meta.data=gex_PS_seurat@meta.data)


gex_combined <- merge(gex_scQer2,y=gex_PS2)



# # # # # # # # # # # # # # # # # # # # 
# # # # #  Integrate and label transfer
# # # # # # # # # # # # # # # # # # # # 



sample_oi <- "Pijuan_Sala_2019"

gex_combined_sub <- subset(gex_combined,
                           subset= ((sample==sample_oi)) | (sample=="mEB_d21_scQer"))


# # # # #  INTEGRATION (following Stuart et al, 2019)
print("splitting combined object")
gex.list <- SplitObject(gex_combined_sub, split.by = "sample")
mem_used()


# normalize and identify variable features for each dataset independently
print("normalizing and finding variable features (separately)")
gex.list <- lapply(X = gex.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
print("selecting integration features")
features_for_integration <- SelectIntegrationFeatures(object.list = gex.list)

print("finding integration anchors")
gex.anchors <- FindIntegrationAnchors(object.list = gex.list, anchor.features = features_for_integration)


print("performing integration")
gex_integrated<- IntegrateData(anchorset = gex.anchors)
DefaultAssay(gex_integrated) <- "integrated"
mem_used()


# Run the standard workflow for visualization and clustering
print("run dim reduction")

n_pcs <- 30
gex_integrated <- ScaleData(gex_integrated, verbose = FALSE)
gex_integrated <- RunPCA(gex_integrated, npcs = n_pcs, verbose = FALSE)
gex_integrated <- RunUMAP(gex_integrated, reduction = "pca", dims = 1:n_pcs)
gex_integrated <- FindNeighbors(gex_integrated, reduction = "pca", dims = 1:n_pcs)
gex_integrated <- FindClusters(gex_integrated, resolution = 0.5)
mem_used()


print("saving integrated out")
saveRDS(gex_integrated,
        sprintf("FindAnchors_integrated_mEB_scQer_vs_%s_20220908.RDS",sample_oi))



# # # # #  LABEL TRANSFER (inspired by Rhodes et al, eLife 2022)
gex_int_sub <- gex_integrated      
PS_ids <-rownames(gex_int_sub@meta.data[gex_int_sub@meta.data$sample == sample_oi,])
cell_oi_ids <- rownames(gex_int_sub@meta.data[gex_int_sub@meta.data$sample == "mEB_d21_scQer",])


print("get integrration embedding")
position_embeds <- gex_int_sub@reductions$pca@cell.embeddings
cells_embeds <- rownames(position_embeds)

pos_PS <- position_embeds[cells_embeds %in% c(PS_ids), ]

maxann_top_nearest <- vector(mode="numeric",length(cell_oi_ids))
names(maxann_top_nearest) <- cell_oi_ids
mostcommon_ann <- vector(mode="character",length(cell_oi_ids))
names(mostcommon_ann) <- cell_oi_ids
list_cell_anno <- list()
list_cell_dist <- list()

n_top <- 100
n_min_top <- 60
print("get nearest neighbours")
for (cell_oi in cell_oi_ids){
  
  print(cell_oi)
  
  pos_oi <- position_embeds[cell_oi, ]
  embed_dist<- sqrt(colSums((t(pos_PS) - pos_oi)^2))
  
  
  # distance to PS from cells in mEB dataset, get closest ones
  embed_dist2 <- embed_dist[order(embed_dist)]
  top_cells<- names(embed_dist2[1:n_top])
  top_embed_dist <- embed_dist2[1:n_top]
  
  #get the annotations of the nearest 5 reference cells
  top_cell_anno <- gex_int_sub@meta.data$celltype[rownames(gex_int_sub@meta.data) %in% top_cells]
  list_cell_anno[[cell_oi]] <- top_cell_anno
  
  list_cell_dist[[cell_oi]] <- top_embed_dist
  
  # mostcommon.ann<- NULL
  # maxann.top_nearest<- NULL
  
  
  #if/else at least 60% match annotations
  maxann<- max(table(top_cell_anno))
  finalann<- names(which.max(table(top_cell_anno)))
  
  maxann_top_nearest[cell_oi]<- maxann
  
  if(maxann >= n_min_top){
    mostcommon_ann[cell_oi]<- finalann
  } else {
    mostcommon_ann[cell_oi]<- "uncertain"
  }
  
}

# generate data frame with cell IDs (consensus annotation, and # top cells with same anno)
df_top_cells <- t(as.data.frame(list_cell_anno))
colnames(df_top_cells) <- paste0("top_anno_",seq(1:n_top))

df_top_cells_dist <- t(as.data.frame(list_cell_dist))
colnames(df_top_cells_dist) <- paste0("top_cell_dist_",seq(1:n_top))

# generate data frame with cell IDs (consensus annotation, and # top cells with same anno)
df_anno <- data.frame(cells=cell_oi_ids,
                      most_common_annotation=mostcommon_ann,
                      number_cells_top_annotation=maxann_top_nearest) %>% 
  cbind(df_top_cells) %>% 
  cbind(df_top_cells_dist)



print("saving annotation label transfer output")  
saveRDS(df_anno,sprintf("/net/shendure/vol1/home/lalannej/vol10_projects_JB/integration_scRNA_mEB_202208/nobackup/annotation_label_transfer_nn_%s_FindAnchors_integration_20220908.RDS",sample_oi))


# }





