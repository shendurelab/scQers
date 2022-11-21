
library(tidyverse)
library(Matrix)
library(stringr)



# load full raw data (including 1s)
df_oBC <- read.table("table_oBC_counts_w_gex_all_reps_v2_all_umis_gex_valid_20220824.txt",
                     header=TRUE)

# gene expression data
GEx_obj_final <-readRDS("GEx_mEB_final_20221025.RDS")
valid_cells <- Cells(GEx_obj_final)
df_oBC <- df_oBC %>% filter(full_cellBC_pdT %in% valid_cells) %>% transform(biol_rep=substr(rep_id,1,str_length(rep_id)-1))


out_dir <- "output_dir/"
setwd("raw_clones_out/")
files <- dir(path="raw_clones_out/",
             pattern=glob2rx("*.RData"))

# parse string
file_names_parsed <- str_split(files,"_")
rep_ids <- lapply(file_names_parsed,"[[",6) %>% unlist

date_str <- "20221027"

# want the unfiltered matrix
oBC_thresh <- 11

# additional thresholds
cells_thresh <- 2
oBC_presence_thresh <- 0.5

rep_ids <- c("A","B","2B")

for (id in seq(length(rep_ids))){
    
  
  rep_oi <-  rep_ids[id]
  print(sprintf("%s",rep_oi))
  
  
  # load raw clonotype data
  load(files[id])
  
  # create output directory
  dir_out_oi <- sprintf("%s/%s/",out_dir,rep_oi)
  out_dir_oBC_ct <- sprintf("%s/%s/oBC_count_submatrices_raw_clones",out_dir,rep_oi)
  
  system(sprintf("mkdir %s",dir_out_oi))
  system(sprintf("mkdir %s",out_dir_oBC_ct))
  
  # get unfiltered oBC sub-matrix
  df_oBC_valid <- df_oBC %>% filter( (biol_rep %in% c(rep_oi)) & UMIs_oBC>=oBC_thresh & p25_valid_oBC & CRE_class=="devCRE")
  
  #valid oBC (list of oBC with at least one cell with above threshold)
  all_valid_oBC <- df_oBC_valid %>% pull(oBC) %>% unique
  all_valid_cells <- df_oBC_valid %>% pull(full_cellBC_pdT) %>% unique()
  
  df_oBC_unfiltered <- df_oBC %>% filter( (full_cellBC_pdT %in% all_valid_cells) & (oBC %in% all_valid_oBC) ) 
  
  cBC_list <- data.frame(cBC_id=seq(length(all_valid_cells)),full_cellBC_pdT=all_valid_cells)
  oBC_list <- data.frame(oBC_id=seq(length(all_valid_oBC)),oBC=all_valid_oBC)
  df_oBC_valid2 <-  df_oBC_unfiltered %>% left_join(cBC_list) %>% left_join(oBC_list)
  
  oBC_mat_unfiltered <- sparseMatrix(j=df_oBC_valid2$oBC_id, i=df_oBC_valid2$cBC_id, x=df_oBC_valid2$UMIs_oBC)
  rownames(oBC_mat_unfiltered) <- cBC_list$full_cellBC_pdT
  colnames(oBC_mat_unfiltered) <- oBC_list$oBC
  oBC_mat_binary <- oBC_mat_unfiltered>=oBC_thresh
  
  
  max_frac_oBC <- c()
  df_clone_info <- data.frame()
  df_oBC_clones <- data.frame()
  
  names_clones <- names(clone_list)
  
  for (name_oi in names_clones){
    
    cells_clone <- clone_list[[name_oi]]
    
    # remove final GEx filtered cells
    cells_clone <- cells_clone[cells_clone %in% valid_cells]
    
    clone_id <- which(name_oi==names_clones)
    
    if (length(cells_clone)>=cells_thresh){
      
      # print(name_oi)
      
      all_oBC_clone <- c()
      df_clone <- data.frame()
      for (cell_oi in cells_clone){
        oBC_list_cell_oi <- colnames(oBC_mat_unfiltered)[which(oBC_mat_binary[cell_oi,])]
        all_oBC_clone <- union(all_oBC_clone,oBC_list_cell_oi)
        
        df_clone <- rbind(df_clone,data.frame(cell=cell_oi,n_oBC=length(oBC_list_cell_oi)))
      }
      
      frac_cell_assigned_w_oBC <- vector(mode="numeric",length=length(all_oBC_clone))
      names(frac_cell_assigned_w_oBC) <- all_oBC_clone
      for (oBC_oi in all_oBC_clone){
        frac_cell_assigned_w_oBC[oBC_oi] <- sum(oBC_mat_binary[cells_clone,oBC_oi])/length(cells_clone)
      }
      
      if (sum(frac_cell_assigned_w_oBC>oBC_presence_thresh)==0){
        next
      } else {

        df_oBC_clones <- rbind(df_oBC_clones,
                               data.frame(clone_id,
                                          clone_name=name_oi,
                                          oBC_pass=all_oBC_clone[frac_cell_assigned_w_oBC>oBC_presence_thresh]))
        
        max_frac_oBC <- c(max_frac_oBC,max(frac_cell_assigned_w_oBC))
        
        df_clone_info <- rbind(df_clone_info,
                               data.frame(clone_id=clone_id,
                                          clone_name=name_oi,
                                          n_cells=length(cells_clone),
                                          n_oBC_25_to_75=sum(frac_cell_assigned_w_oBC>0.25 & frac_cell_assigned_w_oBC<0.75),
                                          n_oBC_75_to_100=sum(frac_cell_assigned_w_oBC>0.75),
                                          max_frac_cell_assigned_w_oBC=max(frac_cell_assigned_w_oBC)))
        # printing raw count output
        mat_oi <- oBC_mat_unfiltered[cells_clone,]
        sum_oBC <- colSums(mat_oi)
        
        idx_sort <- sort(sum_oBC,index.return=TRUE,decreasing=TRUE)
        
        sorted_mat_oi <- t(mat_oi[,idx_sort$ix]) %>% as.matrix
        frac_cells_w_oBC <-  colSums(mat_oi[,idx_sort$ix]>oBC_thresh)/length(cells_clone)
        
        mean_oBC_count <-  colSums(mat_oi[,idx_sort$ix])/length(cells_clone)
        
        sorted_mat_oi_subset <- sorted_mat_oi[mean_oBC_count>1 /length(cells_clone),]
        
        n_cells_w_oBC <-  colSums(sorted_mat_oi_subset>oBC_thresh)
        idx_cells <- sort(n_cells_w_oBC,index.return=TRUE,decreasing=TRUE)
        sorted_mat_oi_subset2 <- sorted_mat_oi_subset[,idx_cells$ix]
        
        write.table(sorted_mat_oi_subset2,
                    sprintf("%s/clone_%d_%s_%s.txt",out_dir_oBC_ct,clone_id,name_oi,date_str),
                    sep="\t",quote=FALSE,col.names = FALSE)
      }
    }
  }
  
  
  summary_plt <- ggplot(df_clone_info)+
    geom_point(aes(x=n_oBC_25_to_75+0.5,y=n_oBC_75_to_100+0.5,color=log10(n_cells)))+
    geom_abline(slope=1)+geom_abline(slope=2,linetype="dashed")#+scale_x_log10()+scale_y_log10()
  
  pdf(sprintf("%s/summary_clone_refinement_%s_%s.pdf",out_dir,rep_oi,date_str),
      width=6,height=5)
  print(summary_plt)
  dev.off()
  

  # valid clones: where thresholding happen
  df_clone_valid_info <- df_clone_info %>% filter(n_oBC_75_to_100>2*n_oBC_25_to_75 & max_frac_cell_assigned_w_oBC>0.9)
  
  # generate oBC clone list: 
  list_oBC_clonotypes <- list()
  for (clone_id_pass in df_clone_valid_info$clone_id){
    list_oBC_clonotypes[[sprintf("clone_%d",clone_id_pass)]] <- df_oBC_clones %>% filter(clone_id==clone_id_pass) %>% pull(oBC_pass)
  }
  
  saveRDS(list_oBC_clonotypes,sprintf("%s/clonotypes_%s_Fisher_%s.RDS",out_dir,rep_oi,date_str))
}




