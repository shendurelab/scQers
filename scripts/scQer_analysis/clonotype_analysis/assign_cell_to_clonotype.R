library(Seurat)
library(dplyr)
#library(patchwork)
#library(sctransform)
#library(ggplot2)
library(tidyverse)
library(scales)
#library(tictoc)
#library(irlba) 
#library(schex) # might need to install
#library(stringdist)
library(Biostrings)
library(stringr)
library(Seurat)
library(Matrix)
library(reshape2)
library(patchwork)
library(igraph)
library(pracma)


# generate the seurat object for the oBC
# filter and reformat

# load full data
df_oBC <- read.table("table_oBC_counts_w_gex_all_reps_v2_all_umis_gex_valid_20220824.txt",
                     header=TRUE)

GEx_obj_final <-readRDS("GEx_mEB_final_20221025.RDS")
valid_cells <- Cells(GEx_obj_final)

df_oBC <- df_oBC %>% filter(full_cellBC_pdT %in% valid_cells) %>% transform(biol_rep=substr(rep_id,1,str_length(rep_id)-1))


# get list of conditions to go through: 
files <- dir(path="raw_clones_out/",
             pattern=glob2rx("*.RData"))

# parse string
file_names_parsed <- str_split(files,"_")
rep_ids <- lapply(file_names_parsed,"[[",6) %>% unlist



date_file <- "20221027"
date_str <- "20221027"

UMI_thresh <- 10 


clonotype_dir <- "refined_clonotypes"
assign_rate <- data.frame()

for (id in seq(length(rep_ids))){
  
  print(id)
  
  rep_oi <-  rep_ids[id]
  
  df_oBC_valid <- df_oBC %>% filter( (biol_rep %in% c(rep_oi))  & p25_valid_oBC & CRE_class=="devCRE")
  
  cBC_list <- df_oBC_valid %>% pull(full_cellBC_pdT) %>% unique()
  cBC_list <- data.frame(cBC_id=seq(length(cBC_list)),full_cellBC_pdT=cBC_list)
  oBC_list <- df_oBC_valid %>% pull(oBC) %>% unique()
  oBC_list <- data.frame(oBC_id=seq(length(oBC_list)),oBC=oBC_list)
  df_oBC_valid2 <-  df_oBC_valid %>% left_join(cBC_list) %>% left_join(oBC_list)
  oBC_mat <- sparseMatrix(j=df_oBC_valid2$oBC_id, i=df_oBC_valid2$cBC_id, x=df_oBC_valid2$UMIs_oBC)
  rownames(oBC_mat) <- cBC_list$full_cellBC_pdT
  colnames(oBC_mat) <- oBC_list$oBC
  
  # transpose oBC count matrix 
  oBC_mat2 <- t(oBC_mat)
  oBC_obj <- CreateSeuratObject(counts = oBC_mat2, min.cells = 0, min.features = 0)

  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # remove nested clusters
  print("removing nested clonotypes")
  
  # oBC_clonotypes <- readRDS(sprintf("%s/clonotypes_%s_%s_Fisher_20220729.RDS",clonotype_dir,rep_ids[id],cell_type_dict[cell_type_oi+1]))
  oBC_clonotypes <- readRDS(sprintf("%s/clonotypes_%s_Fisher_%s.RDS",clonotype_dir,rep_ids[id],date_file))
  
  
  #remove empty clonotypes (too few cells or otherwise no oBC called)
  list_all_oBC <- oBC_clonotypes
  list_all_oBC <- list_all_oBC[!unlist(lapply(list_all_oBC,is_empty))]
  
  # collapse identical clusters
  dup_bool <- duplicated(list_all_oBC)
  list_all_oBC <- list_all_oBC[!dup_bool]
  bool_nested_in <- matrix(data=0,ncol=length(list_all_oBC),nrow=length(list_all_oBC))
  rownames(bool_nested_in) <- names(list_all_oBC)
  colnames(bool_nested_in) <- names(list_all_oBC)
  for (clus1 in names(list_all_oBC)){
    # print(clus1)
    for (clus2 in names(list_all_oBC)){
      if (clus1!=clus2){
        shared_oBC <- intersect(list_all_oBC[[clus1]], list_all_oBC[[clus2]])
        # boolean to find fully nested cluster
        bool_nested_in[clus1,clus2] <- length(shared_oBC) == length(list_all_oBC[[clus1]])
      }
    }
  }
  
  
  nested_clusters <-  names(list_all_oBC)[which(rowSums(bool_nested_in)>0)]
  print("nested clusters:")
  print(nested_clusters)
  list_all_oBC2 <- list_all_oBC
  list_all_oBC2[nested_clusters] <- NULL
  oBC_clonotypes <- list_all_oBC2
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # assign cellBC to clonotypes or classes (doublets, dropout)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  source('perform_clonotype_assignment_v2_20220623.R')
  
  print("assigning cells")
  sample_cell_oi <- sprintf("%s", rep_oi)
  
  run_name <- "clonotype_Fisher_heuristic"
  

  path_out <- sprintf("final_out/%s",
                      sample_cell_oi)
  
  path_out_oBC_counts_mat <- sprintf("final_out/%s/oBC_counts_mats",
                      sample_cell_oi)
  
  system(sprintf("mkdir %s",path_out))
  system(sprintf("mkdir %s",path_out_oBC_counts_mat))
  
  
  df_cBC_metadata_oBC <- perform_clonotype_assignment_v2_20220623(oBC_clonotypes, 
                                                                  oBC_obj, 
                                                                  UMI_thresh,
                                                                  path_out, 
                                                                  sample_cell_oi,
                                                                  run_name)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # final cell classification, print summary figures and final outputs
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # source("/Users/jbl/Documents/UW/Data/sequencing_runs/scMPRA_pilot_data_seq016_seq017_20220202/clonotype_assignment/summary_plot_classification_oBC_20220218.R")
  source("/Users/jbl/Documents/UW/projects/clonotype_mapping_v2_20220622/clonotype_mapping_promoter_series_final_20220729/summary_plot_classification_oBC_v2_20220729.R")
  n_oBC_thresh <- 1.5
  p_to_top_thresh <- 0.975
  recovered_oBC_frac_thresh <- 0.5
  
  df_cBC_metadata_oBC2 <- summary_plot_classification_oBC_v2_20220729(df_cBC_metadata_oBC,
                                                                   oBC_clonotypes,
                                                                   n_oBC_thresh, 
                                                                   p_to_top_thresh, 
                                                                   recovered_oBC_frac_thresh,
                                                                   path_out,
                                                                   run_name,
                                                                   sample_cell_oi)
  
  
  
  df_MOI <- df_oBC_valid %>% group_by(full_cellBC_pdT) %>% summarize(MOI_oBC=sum(UMIs_oBC>UMI_thresh))
  
  df_cBC_metadata_oBC3 <- df_cBC_metadata_oBC2
  df_cBC_metadata_oBC3$category[df_cBC_metadata_oBC3$category=="missed_clonotype" & df_cBC_metadata_oBC3$number_assigned_oBC<=2] <- "missed_clonotype_MOI_less_eq2" 
  df_cBC_metadata_oBC3$category[df_cBC_metadata_oBC3$category=="missed_clonotype" & df_cBC_metadata_oBC3$number_assigned_oBC>2] <- "missed_clonotype_MOI_more2" 
  
  # number singleton clones
  df_cell_counts_per_clone <- df_cBC_metadata_oBC3 %>% filter(category=="singlet") %>% group_by(top_oBC_clone) %>% summarize(n_cells=length(cBC))
  singleton_clones <- df_cell_counts_per_clone %>% filter(n_cells==1) %>% pull(top_oBC_clone)
  
  df_cBC_metadata_oBC4 <- df_cBC_metadata_oBC3
  df_cBC_metadata_oBC4$category[df_cBC_metadata_oBC4$category=="singlet" & (df_cBC_metadata_oBC4$top_oBC_clone %in% singleton_clones)] <- "singlet_singleton"
  
  
  write.table(df_cBC_metadata_oBC4,sprintf("%s/cellBC_assignment_to_clonotypes_%s_v2_clonotype_Fisher_heuristic_%s.txt",path_out,sample_cell_oi,date_str),
              quote=FALSE,sep="\t",row.names=FALSE)
  

  assign_rate <- rbind(assign_rate, data.frame(sample_cell_oi,pivot_wider(as.data.frame(table(df_cBC_metadata_oBC4$category)),names_from=Var1,values_from=Freq)))

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # print count matrices from assignments
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  source("print_oBC_final_clone_assignment_20220729.R")
  
  print_oBC_final_clone_assignment_20220729(oBC_mat,
                                            oBC_clonotypes,
                                            df_cBC_metadata_oBC2,
                                            UMI_thresh,
                                            path_out_oBC_counts_mat,
                                            sample_cell_oi)
    
}




