
# library(Seurat)
library(tidyverse)
library(patchwork)
# library(ArchR)
# library(RColorBrewer)
# library(parallel)
# library(ggrepel)
library(psych)
# library(pryr)



# take in argument
args <- commandArgs(trailingOnly=TRUE)

lookup_str <- args[1]
CRE_oi <- str_split(lookup_str,"__")[[1]][1]
CRE_type <- str_split(lookup_str,"__")[[1]][2]
cluster_group_type <- str_split(lookup_str,"__")[[1]][3]

n_bootstrap <- args[2]
winsor_cut <- args[3]
date_str <- args[4]
path_data <- args[5]


# data: mapped_devCRE_oBC_mBC_mEB_all_reps_v2_UMI_oBC_10_20220822.RDS
df_w_mapped_CRE <- readRDS(path_data) # cell barcode with mapped scQer triplet data (CRE with oBC and mBC UMI counts) 

mean_gex <- df_w_mapped_CRE %>% select(full_cellBC_pdT,norm_gex_UMI) %>% unique() %>% pull(norm_gex_UMI) %>% mean()



options(dplyr.summarise.inform = FALSE)

biol_rep_list <- c("A","B","2B")


# what to average over? 
if (cluster_group_type=="coarse"){
  coarse_clusters <- list(c(0,5),c(1,7),c(4,6),c(2,9,12),3)
  names(coarse_clusters) <- c("pluri","ecto","endo","meso","unk")
} else if (cluster_group_type=="unique"){
  coarse_clusters <- as.list(seq(0,9))
  names(coarse_clusters) <- paste0("unique_",seq(0,9))
} else if (cluster_group_type=="pairs"){
  clusters_included <- seq(0,9)
  coarse_clusters <- combn(clusters_included, 2) %>% as.data.frame %>% as.list
  names(coarse_clusters) <- sprintf("pairs_%d_%d",
                                    coarse_clusters %>% lapply("[[",1) %>% unlist,
                                    coarse_clusters %>% lapply("[[",2) %>% unlist)
}


# assign raw high resolution clusters to group 
cluster_ids <- df_w_mapped_CRE$gex_cluster_id
coarse_cluster_ids <- vector(mode="character",length=length(cluster_ids))
for (coarse_cluster_oi in names(coarse_clusters)){
  cluster_idx <- which(cluster_ids %in% coarse_clusters[[coarse_cluster_oi]])
  coarse_cluster_ids[cluster_idx] <- coarse_cluster_oi
}
df_w_mapped_CRE3 <- df_w_mapped_CRE %>% transform(coarse_cluster_id=coarse_cluster_ids) %>% filter(coarse_cluster_id!="")
df_w_mapped_CRE_subset <- df_w_mapped_CRE3 %>% filter(CRE_id %in% CRE_oi)


df_n_int <- df_w_mapped_CRE3 %>% filter(CRE_class=="devCRE") %>% group_by(biol_rep,CRE_id) %>% summarize(n_integration=n())
df_n_int_median <- df_n_int %>% group_by(biol_rep) %>% summarize(median_int=median(n_integration))


df_CRE_summary_all_boot_CRE_oi <- data.frame()
for (biol_r in biol_rep_list){
  
  df_w_mapped_CRE_subset2 <- df_w_mapped_CRE_subset %>% filter(biol_rep==biol_r)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # # # # NON-PERMUTED
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  print(sprintf("Non permuted: %s, %s", CRE_oi, biol_r))
  
  df_w_mapped_CRE_subset_boot <- df_w_mapped_CRE_subset2 
  df_CRE_summary_by_cluster <- df_w_mapped_CRE_subset_boot %>% group_by(CRE_id,CRE_class,biol_rep,coarse_cluster_id)  %>% 
    summarize(n_integrations=length(unique(cellBC_oBC)),
              n_cells_w_integration=length(unique(full_cellBC_pdT)),
              mean_UMI_per_mBC=mean(UMIs_mBC),
              mean_norm_UMI_per_mBC=mean(UMIs_mBC/norm_gex_UMI*mean_gex))
  
  
  
  major_clusters <- names(coarse_clusters) 
  df_CRE_summary_by_not_cluster <- data.frame()
  for (cluster_oi in major_clusters){
    # print(cluster_oi)
    df_CRE_summary_by_not_cluster_oi <- df_w_mapped_CRE_subset_boot %>% filter(coarse_cluster_id != cluster_oi) %>%
      group_by(CRE_id,biol_rep)  %>% 
      summarize(n_integrations_not_cluster=length(unique(cellBC_oBC)),
                n_cells_w_integration_not_cluster=length(unique(full_cellBC_pdT)),
                mean_UMI_per_mBC_not_cluster=mean(UMIs_mBC),
                mean_norm_UMI_per_mBC_not_cluster=mean(UMIs_mBC/norm_gex_UMI*mean_gex))
    
    df_CRE_summary_by_not_cluster_oi2 <- df_CRE_summary_by_not_cluster_oi %>% transform(coarse_cluster_id=cluster_oi)
    df_CRE_summary_by_not_cluster <- rbind(df_CRE_summary_by_not_cluster,df_CRE_summary_by_not_cluster_oi2)
  }
  
  
  df_CRE_summary_v2 <- df_CRE_summary_by_cluster %>% left_join(df_CRE_summary_by_not_cluster,by = c("CRE_id", "biol_rep", "coarse_cluster_id")) %>%
    transform(FC_cluster_mBC_UMI = mean_UMI_per_mBC/mean_UMI_per_mBC_not_cluster,
              FC_cluster_mBC_norm_UMI = mean_norm_UMI_per_mBC/mean_norm_UMI_per_mBC_not_cluster)
  
  
  
  df_CRE_summary_v3 <- df_CRE_summary_v2 %>% group_by(biol_rep,CRE_id,CRE_class) %>%
    summarize(max_expression_cluster_id = coarse_cluster_id[mean_UMI_per_mBC==max(mean_UMI_per_mBC)][1],
              max_cluster_mBC_UMI = max(mean_UMI_per_mBC),
              max_cluster_FC_mBC_UMI = FC_cluster_mBC_UMI[mean_UMI_per_mBC==max(mean_UMI_per_mBC)][1],
              max_norm_expression_cluster_id = coarse_cluster_id[mean_norm_UMI_per_mBC==max(mean_norm_UMI_per_mBC)][1],
              max_cluster_norm_mBC_UMI = max(mean_norm_UMI_per_mBC),
              max_cluster_norm_FC_mBC_UMI = FC_cluster_mBC_norm_UMI[mean_norm_UMI_per_mBC==max(mean_norm_UMI_per_mBC)][1]) %>% 
    transform(boot_id=NA,permuted=FALSE)
  
  df_CRE_summary_all_boot_CRE_oi <- rbind(df_CRE_summary_all_boot_CRE_oi,df_CRE_summary_v3)
  
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # # # # PERMUTED
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  df_clusters <- df_w_mapped_CRE_subset2 %>% dplyr::select(full_cellBC_pdT,coarse_cluster_id) %>% unique()
  
  
  for (boot in seq(1,n_bootstrap)){
    
    print(sprintf("permuted: %d, %s", boot, biol_r))
    
    df_clusters_reshuf <- df_clusters %>% transform(coarse_cluster_id=sample(coarse_cluster_id,length(coarse_cluster_id), replace=FALSE))
    df_w_mapped_CRE_subset2_reshuf <- df_w_mapped_CRE_subset2 %>% dplyr::select(-coarse_cluster_id) %>% left_join(df_clusters_reshuf)
    
    
    # devCRE refers to elements taken from mouse EB. promoters to exogenous promotesrs
    if (CRE_type=="devCRE"){
      df_w_mapped_CRE_subset_boot <- df_w_mapped_CRE_subset2_reshuf 
      
    } else if (CRE_type=="promoters"){
      n_integration_oi <- df_n_int_median %>% filter(biol_rep==biol_r) %>% pull(median_int) #n_integrations %>% filter(biol_rep==biol_r & CRE_id==CRE_oi) %>% pull(n_int)
      idx_CRE_oi <- which(df_w_mapped_CRE_subset2$CRE_id==CRE_oi)
      sampling_CRE_oi <-  sample(idx_CRE_oi, size=n_integration_oi, replace = TRUE)
      df_w_mapped_CRE_subset_boot <- df_w_mapped_CRE_subset2_reshuf[sampling_CRE_oi,]
    }
    
    df_CRE_summary_by_cluster <- df_w_mapped_CRE_subset_boot %>% group_by(CRE_id,CRE_class,biol_rep,coarse_cluster_id)  %>% 
      summarize(n_integrations=length(unique(cellBC_oBC)),
                n_cells_w_integration=length(unique(full_cellBC_pdT)),
                mean_UMI_per_mBC=mean(UMIs_mBC),
                mean_norm_UMI_per_mBC=mean(UMIs_mBC/norm_gex_UMI*mean_gex))
    
    
    major_clusters <- names(coarse_clusters) 
    df_CRE_summary_by_not_cluster <- data.frame()
    for (cluster_oi in major_clusters){
      # print(cluster_oi)
      df_CRE_summary_by_not_cluster_oi <- df_w_mapped_CRE_subset_boot %>% filter(coarse_cluster_id != cluster_oi) %>%
        group_by(CRE_id,biol_rep)  %>% 
        summarize(n_integrations_not_cluster=length(unique(cellBC_oBC)),
                  n_cells_w_integration_not_cluster=length(unique(full_cellBC_pdT)),
                  mean_UMI_per_mBC_not_cluster=mean(UMIs_mBC),
                  mean_norm_UMI_per_mBC_not_cluster=mean(UMIs_mBC/norm_gex_UMI*mean_gex))
      
      df_CRE_summary_by_not_cluster_oi2 <- df_CRE_summary_by_not_cluster_oi %>% transform(coarse_cluster_id=cluster_oi)
      df_CRE_summary_by_not_cluster <- rbind(df_CRE_summary_by_not_cluster,df_CRE_summary_by_not_cluster_oi2)
    }
    
    
    df_CRE_summary_v2 <- df_CRE_summary_by_cluster %>% left_join(df_CRE_summary_by_not_cluster,by = c("CRE_id", "biol_rep", "coarse_cluster_id")) %>%
      transform(FC_cluster_mBC_UMI = mean_UMI_per_mBC/mean_UMI_per_mBC_not_cluster,
                FC_cluster_mBC_norm_UMI = mean_norm_UMI_per_mBC/mean_norm_UMI_per_mBC_not_cluster)
    
    
    
    df_CRE_summary_v3 <- df_CRE_summary_v2 %>% group_by(biol_rep,CRE_id,CRE_class) %>%
      summarize(max_expression_cluster_id = coarse_cluster_id[mean_UMI_per_mBC==max(mean_UMI_per_mBC)][1],
                max_cluster_mBC_UMI = max(mean_UMI_per_mBC),
                max_cluster_FC_mBC_UMI = FC_cluster_mBC_UMI[mean_UMI_per_mBC==max(mean_UMI_per_mBC)][1],
                max_norm_expression_cluster_id = coarse_cluster_id[mean_norm_UMI_per_mBC==max(mean_norm_UMI_per_mBC)][1],
                max_cluster_norm_mBC_UMI = max(mean_norm_UMI_per_mBC),
                max_cluster_norm_FC_mBC_UMI = FC_cluster_mBC_norm_UMI[mean_norm_UMI_per_mBC==max(mean_norm_UMI_per_mBC)][1]) %>% 
      transform(boot_id=boot,permuted=TRUE)
    
    df_CRE_summary_all_boot_CRE_oi <- rbind(df_CRE_summary_all_boot_CRE_oi,df_CRE_summary_v3)
    
  }
}

write.table(df_CRE_summary_all_boot_CRE_oi,sprintf("permutation_cluster_id_v3_%s_%s_%s_%s.txt",cluster_group_type,CRE_type,CRE_oi,date_str),
            sep="\t", row.names=FALSE, quote=FALSE)








