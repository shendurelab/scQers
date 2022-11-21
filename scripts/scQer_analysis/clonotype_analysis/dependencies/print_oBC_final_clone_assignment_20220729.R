print_oBC_final_clone_assignment_20220729 <- function (oBC_mat,
                                                       oBC_clonotypes,
                                                       df_cBC_metadata_oBC2,
                                                       oBC_thresh,
                                                       out_dir_oBC_ct,
                                                       run_name){
  
  oBC_mat_unfiltered <- oBC_mat
  
  for (clone_oi in names(oBC_clonotypes)){
    
    cells_clonotype <- df_cBC_metadata_oBC2 %>% filter(category=="singlet" & top_oBC_clone==clone_oi) %>% pull(cBC)
    
    if (length(cells_clonotype)>1){
      
      # reps <- str_split(cells_clonotype,"_") %>% lapply("[[",1) %>% unlist
      # cBC <- str_split(cells_clonotype,"_") %>% lapply("[[",2) %>% unlist
      # cells_clone <- paste0(reps,"_",cBC)
      cells_clone <-cells_clonotype
        
      # printing count output
      mat_oi <- oBC_mat_unfiltered[cells_clone,]
      sum_oBC <- colSums(mat_oi)
      
      idx_sort <- sort(sum_oBC,index.return=TRUE,decreasing=TRUE)
      
      sorted_mat_oi <- t(mat_oi[,idx_sort$ix]) %>% as.matrix
      frac_cells_w_oBC <-  colSums(mat_oi[,idx_sort$ix]>oBC_thresh)/length(cells_clone)
      
      mean_oBC_count <-  colSums(mat_oi[,idx_sort$ix])/length(cells_clone)
      
      sorted_mat_oi_subset <- sorted_mat_oi[mean_oBC_count>1/length(cells_clone),]
      
      n_cells_w_oBC <-  colSums(sorted_mat_oi_subset>oBC_thresh)
      idx_cells <- sort(n_cells_w_oBC,index.return=TRUE,decreasing=TRUE)
      sorted_mat_oi_subset2 <- sorted_mat_oi_subset[,idx_cells$ix]
      
      date_str <- str_replace_all(Sys.Date(),"-","") 
      
      write.table(sorted_mat_oi_subset2,
                  sprintf("%s/clone_%s_%s_%s.txt",out_dir_oBC_ct,clone_oi,run_name,date_str),
                  sep="\t",quote=FALSE)
      
    }
  }
}