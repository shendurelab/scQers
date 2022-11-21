
library(tidyverse)
library(stringr)
library(Matrix)


# read oBC data 
df_oBC <- read.table("table_oBC_counts_w_gex_all_reps_umi2_gex_valid_v2_20220823.txt",header=TRUE)
rep_oi <- "A"
permute_id <- 1 
thresh_oBC <- 11


# valid oBC (list of above threshold oBC in the list of )
df_oBC_valid <- df_oBC %>% filter( (biol_rep %in% c(rep_oi)) & UMIs_oBC>=thresh_oBC & p25_valid_oBC & CRE_class=="devCRE") # & gex_cluster_id==cell_type_oi)


# generate oBC count matrix
cBC_list <- df_oBC_valid %>% pull(full_cellBC_pdT) %>% unique()
cBC_list <- data.frame(cBC_id=seq(1:length(cBC_list)),full_cellBC_pdT=cBC_list)
oBC_list <- df_oBC_valid %>% pull(oBC) %>% unique()
oBC_list <- data.frame(oBC_id=seq(1:length(oBC_list)),oBC=oBC_list)

df_oBC_valid2 <-  df_oBC_valid %>% left_join(cBC_list) %>% left_join(oBC_list)

oBC_mat <- sparseMatrix(j=df_oBC_valid2$oBC_id, i=df_oBC_valid2$cBC_id, x=df_oBC_valid2$UMIs_oBC)
rownames(oBC_mat) <- cBC_list$full_cellBC_pdT
colnames(oBC_mat) <- oBC_list$oBC


# threshold for assignment
oBC_mat_binary <- oBC_mat>=thresh_oBC

clone_list <- list()
n_oBC <- dim(oBC_mat)[2]
n_cells <- dim(oBC_mat)[1]
p_val_thresh <- 0.05/(n_cells*n_cells/2)
cell_category <- vector(mode="character",length=n_cells)


# if we want to have multiple repeat of the algorithm with different cell orders, can shuffle order 
shuffled_cell_ids <- sample(rownames(oBC_mat),size=length(rownames(oBC_mat)),replace=FALSE)
names(cell_category) <- shuffled_cell_ids
  
counter <- 1

for (cell_id in shuffled_cell_ids){
  
  if (length(clone_list)==0){
    clone_list[cell_id] <- cell_id
    
    # clone_list <- append(clone_list,cell_id)
    cell_category[cell_id] <- "new_clone"
    
  } else {
    
    
    oBC_list_cell_oi <- which(oBC_mat_binary[cell_id,])
    n_oBC_cell_oi <- length(oBC_list_cell_oi)
    
    member_clones <- c()
    for (clone_oi in names(clone_list)){
      
      oBC_list_clone <- which(oBC_mat_binary[clone_oi,])
      
      n_oBC_clone_oi <- length(oBC_list_clone)
      
      n_in_both <- length(intersect(oBC_list_cell_oi,oBC_list_clone))
      n_in_either <- length(union(oBC_list_cell_oi,oBC_list_clone))
      
      # build contingency table
      n_in_clone_not_in_cell <- n_oBC_clone_oi-n_in_both
      n_in_cell_not_in_clone <- n_oBC_cell_oi-n_in_both
      
      # one-sided Fisher's exact test
      p_val <- fisher.test(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_oBC-n_in_either)), alternative="greater")$p.value
      
      if (p_val<p_val_thresh){
        member_clones <- c(member_clones,clone_oi)
      }
    }
    
    
    if (length(member_clones)==0){
      # no overlap found
      clone_list[cell_id] <- cell_id
      
      cell_category[cell_id] <- "new_clone"
      print(sprintf('%d, cell_id: %s, category: %s',counter, cell_id,"new_clone"))
      
    } else if (length(member_clones)==1){
      # one overlap found: plausible singlet
      clone_list[[member_clones]] <- c(clone_list[[member_clones]],cell_id)
      cell_category[cell_id] <- "singlet"
      
      print(sprintf('%d, cell_id: %s, category: %s',counter,cell_id,"singlet"))

    } else if (length(member_clones)>1){
      cell_category[cell_id] <- "doublet"
      print(sprintf('%d, cell_id: %s, category: %s',counter,cell_id,"doublet"))
      
    }
  }
  counter <- counter+1
}

# save output:
save("clone_list","cell_category","shuffled_cell_ids",
        file=sprintf("clone_list_Fisher_mEBs_rep_%s_permute_id%d.RData",rep_oi,permute_id))

