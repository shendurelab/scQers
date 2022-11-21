perform_clonotype_assignment_v2_20220623 <- function(oBC_clonotypes, 
                                                     oBC_obj2, 
                                                     UMI_thresh,
                                                     path_outs, 
                                                     sample_name, 
                                                     run_name){
  
  
  
  clone_clusters <- names(oBC_clonotypes)
  
  
  # variable initialization
  all_cBCs2 <- Cells(oBC_obj2)
  frac_assigned_oBC_to_clone <- matrix(data=NA,
                                       ncol=length(clone_clusters),
                                       nrow=length(all_cBCs2))
  rownames(frac_assigned_oBC_to_clone) <- all_cBCs2
  colnames(frac_assigned_oBC_to_clone) <- clone_clusters
  
  frac_recovered_oBC_clone <- matrix(data=NA,
                                     ncol=length(clone_clusters),
                                     nrow=length(all_cBCs2))
  rownames(frac_recovered_oBC_clone) <- all_cBCs2
  colnames(frac_recovered_oBC_clone) <- clone_clusters
  
  max_frac_to_cluster <- vector(mode="numeric",length=length(all_cBCs2))
  names(max_frac_to_cluster) <- all_cBCs2
  second_max_frac_to_cluster <- vector(mode="numeric",length=length(all_cBCs2))
  names(second_max_frac_to_cluster) <- all_cBCs2
  number_assigned_oBC <- vector(mode="numeric",length=length(all_cBCs2))
  names(number_assigned_oBC) <- all_cBCs2
  frac_recovered_oBC_top_clone <- vector(mode="numeric",length=length(all_cBCs2))
  names(frac_recovered_oBC_top_clone) <- all_cBCs2
  top_oBC_cluster <- vector(mode="character",length=length(all_cBCs2))
  names(top_oBC_cluster) <- all_cBCs2
  
  
  clusters_included <- as.character(clone_clusters)
  
  sp_BC_count2 <- oBC_obj2@assays$RNA@counts
  count <- 1
  for (cBC_oi in all_cBCs2){
    
    
    sp_count_mat_oi <- sp_BC_count2[,cBC_oi]
    assigned_oBC <- names(which(sp_count_mat_oi>UMI_thresh))
    
    frac_assigned_oBC_to_clone_cBC_oi <- vector(mode="numeric",length=length(clusters_included))
    names(frac_assigned_oBC_to_clone_cBC_oi) <- clusters_included
    
    frac_clone_oBC_from_assigned_cBC_oi <- vector(mode="numeric",length=length(clusters_included))
    names(frac_clone_oBC_from_assigned_cBC_oi) <- clusters_included
    
    for (cluster_oi in clusters_included){
      
      # fraction of eBC assigned to cell to core eBC of high confidence clonotypes
      frac_assigned_oBC_to_clone_cBC_oi[cluster_oi] <- length(intersect(assigned_oBC,oBC_clonotypes[[cluster_oi]]))/length(assigned_oBC)
      frac_clone_oBC_from_assigned_cBC_oi[cluster_oi] <- length(intersect(assigned_oBC,oBC_clonotypes[[cluster_oi]]))/length(oBC_clonotypes[[cluster_oi]])
      
    }
    number_assigned_oBC[cBC_oi] <- length(assigned_oBC)
    
    if (number_assigned_oBC[cBC_oi]>0){
      
      sorted_frac_core_oBC_cluster <- sort(frac_assigned_oBC_to_clone_cBC_oi,decreasing=TRUE)
      
      max_frac_to_cluster[cBC_oi] <- sorted_frac_core_oBC_cluster[1]
      frac_recovered_oBC_top_clone[cBC_oi] <- frac_clone_oBC_from_assigned_cBC_oi[names(sorted_frac_core_oBC_cluster[1])]
      
      if (max_frac_to_cluster[cBC_oi]>0){
        top_clono <- names(which(frac_assigned_oBC_to_clone_cBC_oi==max(frac_assigned_oBC_to_clone_cBC_oi)))
        if (length(top_clono)==1){
          top_oBC_cluster[cBC_oi] <- top_clono
        }
        else{
          top_oBC_cluster[cBC_oi] <- top_clono[1]
          # print(sprintf("tie at the top for %s",cBC_oi))
        }
      }
    
      second_max_frac_to_cluster[cBC_oi] <- sorted_frac_core_oBC_cluster[2]
    }
    
    
    count <- count+1
    if (( count %% 300)==0){
      print(count)
    }
    
    frac_assigned_oBC_to_clone[cBC_oi,clusters_included] <- frac_assigned_oBC_to_clone_cBC_oi[clusters_included]
    frac_recovered_oBC_clone[cBC_oi,clusters_included] <- frac_clone_oBC_from_assigned_cBC_oi[clusters_included]
  }
  
  
  df_cBC_metadata_oBC <- data.frame(cBC=all_cBCs2,
                                    number_assigned_oBC,
                                    frac_oBC_assigned_to_top_clone=max_frac_to_cluster,
                                    frac_oBC_recovered_top_clone=frac_recovered_oBC_top_clone,
                                    top_oBC_clone=top_oBC_cluster,
                                    # original_seurat_cluster=oBC_obj2$seurat_clusters[all_cBCs2],
                                    second_max_frac_to_clone=second_max_frac_to_cluster)
  
  date_str <- str_replace_all(Sys.Date(),"-","") 
  write.table(df_cBC_metadata_oBC,sprintf("%s/cellBC_assignment_to_clonotypes_%s_%s_%s.txt",path_outs,sample_name,run_name,date_str),
              sep="\t",row.names=FALSE,quote=FALSE)
  
  return(df_cBC_metadata_oBC)
  
}