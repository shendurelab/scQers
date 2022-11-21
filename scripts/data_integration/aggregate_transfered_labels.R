

library(Seurat)
library(tidyverse)


# generate coarse-grained cell type table
int_oi <- "FindAnchors"
file_anno <- sprintf("annotation_label_transfer_nn_Pijuan_Sala_2019_%s_integration_20220908.RDS",int_oi)
df_anno <- readRDS(file_anno)
all_cell_types <- c(df_anno$top_anno_1, df_anno$top_anno_2)
all_cell_types2 <- unique(all_cell_types)
df_coarse_group <- data.frame(nn_annotation=all_cell_types2[all_cell_types2 !="uncertain"],
                              nn_coarse_annotation=all_cell_types2[all_cell_types2 !="uncertain"])
primitive_streak_group <- c("Primitive Streak","Anterior Primitive Streak")
endoderm_group <- c("ExE endoderm","Parietal endoderm","Visceral endoderm","gut","Def. endoderm","Gut")
neuroectoderm_group <- c("Caudal neurectoderm","Spinal cord","Forebrain/Midbrain/Hindbrain","Rostral neurectoderm")
mesoderm_group <- c("Caudal epiblast","Nascent mesoderm","Mixed mesoderm","Mesenchyme",
                    "Intermediate mesoderm","Caudal Mesoderm","Paraxial mesoderm",
                    "Pharyngeal mesoderm","Somitic mesoderm","ExE mesoderm")
haemaendothelium_group <- c("Erythroid1","Erythroid3","Erythroid2","Blood progenitors 1","Blood progenitors 2","Endothelium","Haematoendothelial progenitors")
df_coarse_group$nn_coarse_annotation[df_coarse_group$nn_annotation %in% primitive_streak_group] <- "coarse_primitive_streak"
df_coarse_group$nn_coarse_annotation[df_coarse_group$nn_annotation %in% endoderm_group] <- "coarse_endoderm"
df_coarse_group$nn_coarse_annotation[df_coarse_group$nn_annotation %in% mesoderm_group] <- "coarse_mesoderm"
df_coarse_group$nn_coarse_annotation[df_coarse_group$nn_annotation %in% neuroectoderm_group] <- "coarse_neuroectoderm"
df_coarse_group$nn_coarse_annotation[df_coarse_group$nn_annotation %in% haemaendothelium_group] <- "coarse_haemaendothelium"
rownames(df_coarse_group) <- df_coarse_group$nn_annotation



# load GEx data
gex_obj <- read.RDS("GEx_obj_sc_rep_mEB_series.RDS")

all_cell_types <- c()

type_int <- c("FindAnchors")

for (int_oi in type_int){
  print(int_oi)
  
  #  transfered labels table
  file_anno <- sprintf("annotation_label_transfer_nn_Pijuan_Sala_2019_%s_integration_20220908.RDS",int_oi)
  
  
  df_anno <- readRDS(file_anno)
  
  df_anno_long <- df_anno %>% select(cells,most_common_annotation,number_cells_top_annotation,starts_with("top_anno_")) %>% pivot_longer(cols=starts_with("top_anno_")) %>% transform(nn_rank=str_replace(name,"top_anno_","") %>% as.numeric)
  df_anno_long2 <- df_anno_long %>% select(cells,most_common_annotation,number_cells_top_annotation,nn_annotation=value,nn_rank)
  
  
  df_dist_long <- df_anno %>% select(cells,starts_with("top_cell_dist_")) %>% pivot_longer(cols=starts_with("top_cell_dist_")) %>% transform(nn_rank=str_replace(name,"top_cell_dist_","") %>% as.numeric)
  df_dist_long2 <- df_dist_long %>% select(cells,nn_rank,distance_embed=value)
  
  df_anno_long3 <- df_anno_long2 %>% left_join(df_dist_long2)
  
  df_anno_w_coarse <- df_anno_long3 %>% left_join(df_coarse_group) 
  
  # number of nearest neighbors used for final label transfer
  n_top <- 10
  print(n_top)
  n_min_top <- n_top*0.6
  df_common_coarse <- df_anno_w_coarse %>% filter(nn_rank<=n_top) %>% group_by(cells) %>% 
    summarize(number_cells_top_annotation = max(table(nn_annotation)),
              top_annotation = names(which.max(table(nn_annotation))),
              mean_dist_top_nn = mean(distance_embed),
              number_cells_top_coarse_annotation = max(table(nn_coarse_annotation)),
              top_coarse_annotation = names(which.max(table(nn_coarse_annotation))) )
  
  df_common_coarse2 <- df_common_coarse %>% data.frame()
  df_common_coarse2$top_annotation[df_common_coarse2$number_cells_top_annotation<n_min_top] <- "uncertain"
  df_common_coarse2$top_coarse_annotation[df_common_coarse2$number_cells_top_coarse_annotation<n_min_top] <- "uncertain"
  
  rownames(df_common_coarse2) <- df_common_coarse2$cells
  
  
  
  gex_obj3 <- AddMetaData(gex_obj,df_common_coarse2)
  
  
  df_metadata <- gex_obj3@meta.data
  df_cells_by_cluster <- df_metadata %>% group_by(seurat_clusters) %>% summarize(n_cells=length(cells))
  
  
  df_anno_summary <- df_metadata %>% 
    group_by(seurat_clusters,top_annotation) %>% 
    summarize(n_per_anno=length(cells)) %>% 
    left_join(df_cells_by_cluster) %>%
    transform(frac_cells_w_anno=n_per_anno/n_cells) %>% 
    select(seurat_clusters,top_annotation,n_per_anno,frac_cells_w_anno)
  
  df_anno_summary_coarse <- df_metadata %>% 
    group_by(seurat_clusters,top_coarse_annotation) %>% 
    summarize(n_per_anno=length(cells)) %>% 
    left_join(df_cells_by_cluster) %>%
    transform(frac_cells_w_anno=n_per_anno/n_cells) %>% 
    select(seurat_clusters,top_coarse_annotation,n_per_anno,frac_cells_w_anno)
  
  

  # order based on coarse annotation
  df_anno2 <- df_anno_summary %>% filter(!(seurat_clusters %in% c(10,12)))
  
  df_wide1 <- df_anno2 %>% select(-n_per_anno) %>% pivot_wider(values_from=frac_cells_w_anno,names_from=top_annotation)
  df_wide1[is.na(df_wide1)] <- 0
  mat1 <- df_wide1[,2:dim(df_wide1)[2]] %>% as.matrix
  rownames(mat1) <- df_wide1$seurat_clusters
  
  cluster_obj <- hclust(dist(mat1), method = "complete", members = NULL)
  cluster_obj2 <- hclust(dist(t(mat1)), method = "complete", members = NULL)
  
  
  df_order2 <- read.table("fine_cluster_PS_ordering_v2_20220927.txt",
                          header=TRUE,sep="\t")
  
  df_order2_coarse <- read.table("coarse_cluster_PS_ordering_v2_20220916.txt",
                          header=TRUE,sep="\t")

  
  df_order <- data.frame(order_clust=seq(1,10),
                         seurat_clusters=c("0","5","1","7","3","6","4","2","9","8"))
  
  df_anno2_wd <- df_anno2 %>% select(-n_per_anno) %>% pivot_wider(values_from=frac_cells_w_anno,names_from=top_annotation)
  df_anno2_wd2 <- df_anno2_wd
  df_anno2_wd2[is.na(df_anno2_wd2)] <- 0
  df_anno2_2 <- df_anno2_wd2 %>% pivot_longer(cols=!(seurat_clusters),names_to="top_annotation",values_to="frac_cells_w_anno")
  
  df_anno3 <- df_anno2_2 %>% left_join(df_order) %>% left_join(df_order2)
  df_anno4 <- df_anno3 %>% mutate(seurat_clusters = fct_reorder(seurat_clusters, order_clust),
                                  top_annotation = fct_reorder(top_annotation,order_clust2))
  
  
  # exclude cell-types with little cells transfered 
  rep_thresh <- 0.05
  lowly_rep_cell_type <- df_anno4 %>% group_by(top_annotation) %>% summarize(max_frac_match=max(frac_cells_w_anno)) %>% filter(max_frac_match<rep_thresh) %>% pull(top_annotation)

  plt_labels_fine <- ggplot(df_anno4 %>% filter(!(top_annotation %in% lowly_rep_cell_type))) +
    geom_tile(aes(y=seurat_clusters,x=top_annotation,fill=frac_cells_w_anno,color=frac_cells_w_anno))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme_cowplot()+
    theme(panel.grid.major=element_line(size=0.1,color="grey"),
          axis.text=element_text(size=11),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.title=element_text(size=10))+coord_fixed()
  

  
  df_anno2_coarse_wd <- df_anno_summary_coarse %>% select(-n_per_anno) %>% pivot_wider(values_from=frac_cells_w_anno,names_from=top_coarse_annotation)
  df_anno2_coarse_wd2 <- df_anno2_coarse_wd
  df_anno2_coarse_wd2[is.na(df_anno2_coarse_wd2)] <- 0
  df_anno2_coarse_2 <- df_anno2_coarse_wd2 %>% pivot_longer(cols=!(seurat_clusters),names_to="top_coarse_annotation",values_to="frac_cells_w_anno")
  
  df_anno3_coarse <- df_anno2_coarse_2 %>% left_join(df_order) %>% left_join(df_order2_coarse)
  df_anno4_coarse <- df_anno3_coarse %>% mutate(seurat_clusters = fct_reorder(seurat_clusters, order_clust),
                                  top_coarse_annotation = fct_reorder(top_coarse_annotation,order_coarse_clust))
  
  
  # exclude cell-types with little cells transfered 
  rep_thresh <- 0.05
  lowly_rep_cell_type_coarse <- df_anno4_coarse %>% group_by(top_coarse_annotation) %>% summarize(max_frac_match=max(frac_cells_w_anno)) %>% filter(max_frac_match<rep_thresh) %>% pull(top_coarse_annotation)
  
  
  plt_labels_coarse <- ggplot(df_anno4_coarse %>% filter(!(seurat_clusters %in% c(10,12)) & !(top_coarse_annotation %in% lowly_rep_cell_type_coarse))) +
    geom_tile(aes(y=seurat_clusters,x=top_coarse_annotation,fill=frac_cells_w_anno,color=frac_cells_w_anno))+
    theme_cowplot()+
    theme(panel.grid.major=element_line(size=0.1,color="grey"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=11),
          axis.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.title=element_text(size=10))+coord_fixed()
  
  plt_design <- "AAAAAA#BBB"
  pdf(sprintf("matrix_labels_by_cluster_%s_n_top_%d_20220927.pdf",int_oi,n_top),
      height=7,width=16)
  print(plt_labels_fine+plt_labels_coarse+plot_layout(design=plt_design))
  dev.off()
  
  pdf(sprintf("matrix_labels_by_cluster_fine_%s_n_top_%d_20220927.pdf",int_oi,n_top),
      height=5,width=20)
  print(plt_labels_fine)
  dev.off()
  
}


