

# get list of conditions to go through: 
files <- dir(path="raw_clones_out/",
             pattern=glob2rx("*.RData"))

# parse string
file_names_parsed <- str_split(files,"_")
rep_ids <- lapply(file_names_parsed,"[[",6) %>% unlist


clonotype_dir <- "refined_clonotypes/"
assign_rate <- data.frame()

dir_data <- "final_out"
date_str <- "20221027"

n_cells_thresh <- 3
df_cell_metadata_all <- data.frame()
df_clone_metadata_all <- data.frame()

list_all_clones <- list()
id <- 1

for (id in seq(length(rep_ids))){
  
  rep_oi <-  rep_ids[id]
  
  # cell_type_oi <- cell_types[[id]]
  sample_cell_oi <- sprintf("%s", rep_oi)
  print(sample_cell_oi)
  
  file_metadata <- sprintf("%s/%s/cellBC_assignment_to_clonotypes_%s_v2_clonotype_Fisher_heuristic_%s.txt",
                           dir_data,sample_cell_oi,sample_cell_oi,date_str) 
  df_metadata <- read.table(file_metadata,header=TRUE,sep="\t")
  
  df_cell_assignment <- df_metadata %>% filter(category=="singlet") %>% select(Cells=cBC,top_oBC_clone,category)
  valid_clones <- df_cell_assignment %>% pull(top_oBC_clone) %>% unique
  rownames(df_cell_assignment) <- df_cell_assignment$Cells
  
  df_metadata2 <- df_metadata %>% transform(rep_id=rep_oi, top_oBC_clone=paste0(sample_cell_oi,"_",top_oBC_clone))
  df_cell_metadata_all <- rbind(df_cell_metadata_all,df_metadata2)
  
  
  # restrict to cells and oBC in final clonotypes
  file_assignments <- sprintf("%s/%s/oBC_clonotypes_metadata_%s_clonotype_Fisher_heuristic_%s.txt",
                              dir_data,sample_cell_oi,sample_cell_oi,date_str)
  
  df_assignment <- read.table(file_assignments,header=TRUE,sep="\t")
  df_assignment2 <- df_assignment %>% filter(clonotype_id %in% valid_clones & n_assigned_cellBC>=n_cells_thresh)
  
  for (clone_oi in (df_assignment2$clonotype_id)){
    list_all_clones[[sprintf("%s_%s",sample_cell_oi,clone_oi)]] <- df_assignment2 %>% filter(clonotype_id==clone_oi) %>% pull(list_oBC) %>% str_split(",") %>% unlist
  }
  
  df_clone_metadata_all <- rbind(df_clone_metadata_all, 
                                 df_assignment2 %>% mutate(clonotype_id=paste0(sample_cell_oi,"_",clonotype_id)) %>% transform(biol_rep=rep_oi))
  
}

setwd(dir_data)

saveRDS(list_all_clones,"list_final_oBC_clonotypes_all_reps.RDS")
write.table(df_cell_metadata_all,"metadata_cell_final_clonotype_assignments_all_reps.txt",
            sep="\t",quote=FALSE,row.names=FALSE)
write.table(df_clone_metadata_all,"metadata_clones_final_assignments_all_reps_mEB.txt",
            sep="\t",quote=FALSE,row.names=FALSE)



