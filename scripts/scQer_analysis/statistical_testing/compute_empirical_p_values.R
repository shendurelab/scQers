

library(tidyverse)

path_boot_out <- "results_bootstrap"
path_perm_out <- "results_permutation"



biol_reps <- c("A","B","2B")
group_types <-  c("coarse","unique","pairs")

for (group_oi in group_types){
  
  
  # bootstrap and comparison to minP and no promoter
  print("loading data")
  df_boot <- readRDS(sprintf("%s/table_bootstrap_minP_noP_all_CREs_%s_all_reps_v3_20220828.RDS",path_boot_out,group_oi))
  
  
  CRE_list <- df_boot %>% pull(CRE_id) %>% unique 
  CRE_list <- CRE_list[!is.na(CRE_list)]
  
  df_activity_p_value <- data.frame()
  
  for (CRE_oi in CRE_list){
    
    for (biol_r in biol_reps){
      
      print(sprintf("%s, %s",CRE_oi,biol_r))
      df_boot_oi <- df_boot %>% filter(bootstrap_for_CRE==CRE_oi & biol_rep==biol_r)
      
      max_cluster_exp <- df_boot_oi %>% filter(CRE_id==CRE_oi) %>% pull(max_cluster_norm_mBC_UMI)
      
      max_cluster_exp_min_noP <- df_boot_oi %>% filter(CRE_id %in% c("minP","noP")) %>% pull(max_cluster_norm_mBC_UMI)
      
      
      # obtaining information about the most frequent max cluster
      max_exp_cluster_id <- df_boot_oi %>% filter(CRE_id==CRE_oi) %>% pull(max_norm_expression_cluster_id)
      freq_max_cluster_id <- table(max_exp_cluster_id)
      most_freq_max_clus <- names(freq_max_cluster_id)[freq_max_cluster_id==max(freq_max_cluster_id)][1]
      freq_most_freq_max_clus <- mean(max_exp_cluster_id==most_freq_max_clus)
      
      empirical_p_val_oi <- vector(mode="numeric",length=length(max_cluster_exp))
      for (idx in seq(length(max_cluster_exp))){
        empirical_p_val_oi[idx] <- (sum(max_cluster_exp_min_noP>=max_cluster_exp[idx])+1)/(length(max_cluster_exp_min_noP)+1)
      }
      df_activity_p_value_oi <- data.frame(CRE_id=CRE_oi,
                                           biol_rep=biol_r,
                                           q25_bootstrap_max_exp_norm=quantile(max_cluster_exp,0.25),
                                           q50_bootstrap_max_exp_norm=quantile(max_cluster_exp,0.5),
                                           q75_bootstrap_max_exp_norm=quantile(max_cluster_exp,0.75),
                                           empirical_bootstrap_p_val_vs_noP_minP=mean(empirical_p_val_oi),
                                           most_frequent_max_exp_cluster = most_freq_max_clus,
                                           frequency_most_frequent_max_exp_cluster = freq_most_freq_max_clus,
                                           group_type=group_oi)
      df_activity_p_value <- rbind(df_activity_p_value,df_activity_p_value_oi)
    }
    
  }
  
  
  # BH correction
  df_activity_p_value2 <- data.frame()
  for (rep_oi in biol_reps){
    df_activity_p_value2_oi <- df_activity_p_value %>% filter(biol_rep==rep_oi)
    df_activity_p_value2_oi$empirical_bootstrap_p_val_vs_noP_minP_BH_corr <- p.adjust(df_activity_p_value2_oi$empirical_bootstrap_p_val_vs_noP_minP,method="BH")
    
    df_activity_p_value2 <- rbind(df_activity_p_value2,df_activity_p_value2_oi)
  }
  

  write.table(df_activity_p_value2,sprintf("p_val_activity_v3_%s_20220829.txt",group_oi),
              sep="\t",row.names=FALSE,quote=FALSE)
  
}



# # # # # # # # # # # # # # # # # # # # # # 
# # # cell-type specific activity
# # # # # # # # # # # # # # # # # # # # # # 


for (group_oi in group_types){
  
  print(group_oi)
  df_perm <- readRDS(sprintf("%s/table_gex_cluster_permutations_all_CREs_%s_all_reps_v3_20220828.RDS",path_perm_out,group_oi))
  df_boot <- readRDS(sprintf("%s/table_bootstrap_minP_noP_all_CREs_%s_all_reps_v3_20220828.RDS",path_boot_out,group_oi))
  
  df_cell_type_specificity_p_value <- data.frame()
  
  
  for (CRE_oi in CRE_list){
      for (biol_r in biol_reps){
      
      print(sprintf("%s, %s",CRE_oi, biol_r))
      
      df_boot_oi <- df_boot %>% filter(bootstrap_for_CRE==CRE_oi & biol_rep==biol_r)
      df_perm_oi <- df_perm %>% filter(CRE_id==CRE_oi & permuted & biol_rep==biol_r)
      
      df_perm_2 <- df_perm_oi %>% select(CRE_id,max_cluster_norm_mBC_UMI,max_cluster_norm_FC_mBC_UMI,biol_rep,permuted) %>%
        rbind(df_boot_oi %>% filter(CRE_id==CRE_oi) %>% select(CRE_id,max_cluster_norm_mBC_UMI,max_cluster_norm_FC_mBC_UMI,biol_rep) %>% transform(permuted=FALSE))
      
      max_cluster_norm_FC_mBC_UMI_non_perm <- df_perm_2 %>% filter(!permuted) %>% pull(max_cluster_norm_FC_mBC_UMI)
      max_cluster_norm_FC_mBC_UMI_perm <- df_perm_2 %>% filter(permuted) %>% pull(max_cluster_norm_FC_mBC_UMI)
      
      max_cluster_norm_FC_mBC_UMI_non_perm[is.na(max_cluster_norm_FC_mBC_UMI_non_perm)] <- 1
      max_cluster_norm_FC_mBC_UMI_perm[is.na(max_cluster_norm_FC_mBC_UMI_perm)] <- 1
      
      empirical_p_val_oi <- vector(mode="numeric",length=length(max_cluster_norm_FC_mBC_UMI_non_perm))
      for (idx in seq(length(max_cluster_norm_FC_mBC_UMI_non_perm))){
        empirical_p_val_oi[idx] <- (sum(max_cluster_norm_FC_mBC_UMI_perm>=max_cluster_norm_FC_mBC_UMI_non_perm[idx])+1)/(length(max_cluster_norm_FC_mBC_UMI_perm)+1)
        if (is.na(empirical_p_val_oi[idx])){
          print(sprintf("problem with %s %s, id=%d",group_oi,CRE_oi,idx))
        }
      }
      
      
      df_cell_type_specificity_p_value_oi <- data.frame(CRE_id=CRE_oi,
                                                        biol_rep=biol_r,
                                                        q25_bootstrap_FC_norm=quantile(max_cluster_norm_FC_mBC_UMI_non_perm,0.25),
                                                        q50_bootstrap_FC_norm=quantile(max_cluster_norm_FC_mBC_UMI_non_perm,0.5),
                                                        q75_bootstrap_FC_norm=quantile(max_cluster_norm_FC_mBC_UMI_non_perm,0.75),
                                                        empirical_permute_gex_cluster_p_val=mean(empirical_p_val_oi),
                                                        group_type=group_oi)
      
      df_cell_type_specificity_p_value <- rbind(df_cell_type_specificity_p_value,df_cell_type_specificity_p_value_oi)
      
      
      }
    
  }
  
  # BH correction
  df_cell_type_specificity_p_value2 <- data.frame()
  for (rep_oi in biol_reps){
    df_cell_type_specificity_p_value2_oi <- df_cell_type_specificity_p_value %>% filter(biol_rep==rep_oi)
    df_cell_type_specificity_p_value2_oi$empirical_permute_gex_cluster_p_val_BH_corr <- p.adjust(df_activity_p_value2_oi$empirical_permute_gex_cluster_p_val,method="BH")
    
    df_cell_type_specificity_p_value2 <- rbind(df_cell_type_specificity_p_value2,df_cell_type_specificity_p_value2_oi)
  }
  
  
  write.table(df_cell_type_specificity_p_value2,sprintf("p_val_specificity_v3_%s_20220829.txt",group_oi),
              sep="\t", quote=FALSE, row.names=FALSE)
  
}



