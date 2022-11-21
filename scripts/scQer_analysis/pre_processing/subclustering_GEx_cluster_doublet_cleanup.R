
library(Seurat)
library(tidyverse)
library(patchwork)
library(ArchR)
library(RColorBrewer)
library(ggplotify) # might need to install






gex_obj <- readRDS("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_scREA_mEB_bottlenecked_all_reps_A_B_2B_20220610.RDS")


gex_meta <- gex_obj@meta.data
biol_rep <- gex_obj$orig.ident  %>% substr(1,str_length(gex_obj$orig.ident)-1)
gex_obj$biol_rep <- biol_rep

plt1 <- DimPlot(gex_obj, reduction = "umap", label=TRUE,split.by = "biol_rep")
plt2 <- FeaturePlot(object = gex_obj,  features = 'doublet_score',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
plt3 <- ggplot(gex_meta)+geom_boxplot(aes(x=seurat_clusters,group=seurat_clusters,y=doublet_score))
plt_design <- "AAAA
BBCC"
plt1+plt2+plt3 +plot_layout(design=plt_design)


cluster_oi <- 8
cells_cluster_oi <- Cells(gex_obj)[gex_obj$seurat_clusters==cluster_oi]

markers_oi <- FindMarkers(gex_obj,
                          ident.1 =cells_cluster_oi,
                          logfc.threshold=0.5)

markers_oi %>% arrange(desc(avg_log2FC))

cells_cluster_oi2 <- Cells(gex_obj)[gex_obj$seurat_clusters==9]

markers_oi <- FindMarkers(gex_obj,
                          ident.2 = cells_cluster_oi, 
                          ident.1 =cells_cluster_oi2,
            logfc.threshold=0.5)


setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/subclustering_analysis_20220819")


clusters <- c(0,1,2,3,4,5,6,7,8,9,10,12)


df_clusters <- c()

# each loop sub-cluster an original cluster
for (cluster_oi in clusters){
  
  print(cluster_oi)
  
  
  # get only cells in a given cluster
  gex_obj_clus_oi <- subset(gex_obj,subset= seurat_clusters==cluster_oi)
  
  gex_obj_clus_oi <- NormalizeData(gex_obj_clus_oi, normalization.method = "LogNormalize", scale.factor = 10000)
  gex_obj_clus_oi <- FindVariableFeatures(gex_obj_clus_oi, selection.method = "vst", nfeatures = 1000, verbose = TRUE)
  all.genes <- rownames(gex_obj_clus_oi)
  gex_obj_clus_oi <- ScaleData(gex_obj_clus_oi, features = all.genes)
  gex_obj_clus_oi <- RunPCA(gex_obj_clus_oi, features = VariableFeatures(object = gex_obj_clus_oi), verbose = FALSE, npcs = 100)
  
  
  gex_obj_clus_oi_2 <- gex_obj_clus_oi
  top_pc <- 50 
  gex_obj_clus_oi_2 <- FindNeighbors(gex_obj_clus_oi_2, dims = 1:top_pc)
  gex_obj_clus_oi_2 <- FindClusters(gex_obj_clus_oi_2, resolution = 0.5) # before was 0.5
  gex_obj_clus_oi_2 <- RunUMAP(gex_obj_clus_oi_2, dims = 1:top_pc, n.neighbors = 50, seed.use = 42)
  
  # biol_rep <- gex_obj_clus_oi_2$orig.ident  %>% substr(1,str_length(gex_obj_clus_oi_2$orig.ident)-1)
  # gex_obj_clus_oi_2$biol_rep <- biol_rep
  
  gex_meta_sub <- gex_obj_clus_oi_2@meta.data
  plt3 <- ggplot(gex_meta_sub)+geom_boxplot(aes(x=RNA_snn_res.0.5,group=RNA_snn_res.0.5,y=doublet_score))
  
  
  plt1 <- DimPlot(gex_obj_clus_oi_2, reduction = "umap", label=TRUE,split.by = "biol_rep")+labs(title=sprintf("cluster %d",cluster_oi))
  
  # plt1 <- DimPlot(gex_obj_clus_oi_2, reduction = "umap", label=TRUE)
  # FeaturePlot(object = gex_obj_clus_oi_2,  features = 'nCount_RNA',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
  plt2 <- FeaturePlot(object = gex_obj_clus_oi_2,  features = 'doublet_score',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
  
  plt_design <- "AAAA
                 BBCC"
  final_plt <- plt1+plt2+plt3 +plot_layout(design=plt_design)
  
  pdf(sprintf("subclustering_%d_gex_mEB_all_reps_20220819.pdf",cluster_oi),width=9,height=7)
  print(final_plt)
  dev.off()
  
  df_clusters_oi <- data.frame(cells= Cells(gex_obj_clus_oi_2),
                               original_cluster=cluster_oi,
                               subclustering_id=gex_obj_clus_oi_2$RNA_snn_res.0.5)
  
  df_clusters <- rbind(df_clusters,df_clusters_oi)
  
}



# cells to remove from sub clustering (inspection of QC figures generated above)

cells_flagged_ori <- Cells(gex_obj)[gex_obj$seurat_clusters==11]

ori_clus <- c(0,3,3,4,5,6,7,9)
sub_clus <- c(4,6,7,5,5,5,3,3)

cell_flagged_subclust <- c()
for (idx in seq(length(ori_clus))){
  cbc_flag <- df_clusters %>% filter(original_cluster==ori_clus[idx] &
                           subclustering_id==sub_clus[idx]) %>% pull(cells)
  
  cell_flagged_subclust <- c(cell_flagged_subclust,cbc_flag)
  
}

cells_flagged <- c(cells_flagged_ori,cell_flagged_subclust)
cells <- Cells(gex_obj)
cells_no_doub <- cells[!(cells %in% cells_flagged)]

gex_obj_v2<- subset(gex_obj,cells= cells_no_doub)

saveRDS(gex_obj_v2,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_scREA_mEB_bottlenecked_all_reps_A_B_2B_subclus_doub_rem_20220819.RDS")


