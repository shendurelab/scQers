# library(Seurat)
# library(scales)
# library(ggrepel)

# setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/statistical_assessment_CREs_20220630/")

# list_prom <- c("pgk1P","minP","noP","ubcP","eef1aP")



# # # 
# regenerate the combined table of output: 

type_groups <- c("*_coarse_*","*_unique_*","*_pairs_*")
names(type_groups) <-  c("coarse","unique","pairs")


path_out_perm <- "results_permutation"
setwd(path_out_perm)


print("Combining the permutation outs:")

for (group_oi in names(type_groups)){
  
  print(group_oi)
  
  files_perm <- dir(pattern = glob2rx(type_groups[group_oi]))
  
  df_oi <- read.table(files_perm[1],header=TRUE)
  n_rows <- dim(df_oi)[1]
  
  df_perm <- data.frame(matrix(NA, nrow = length(files_perm)*n_rows, ncol = dim(df_oi)[2]))
  colnames(df_perm) <- colnames(df_oi)
  list_problem_CRE <- c()
  counter <- 1
  for (file_oi in files_perm){
    print(file_oi)
    df_oi <- read.table(file_oi,header=TRUE)
    idx <- seq(n_rows)+(counter-1)*n_rows
    if (dim(df_oi)[1]!=n_rows){
      print(sprintf("problem with %s",file_oi))
      list_problem_CRE <- c(list_problem_CRE,file_oi)
    } else {
      df_perm[idx,] <- df_oi
    }
    counter <- counter+1
  }
  
  saveRDS(df_perm,sprintf("table_gex_cluster_permutations_all_CREs_%s_all_reps_v3_20220823.RDS",group_oi))
  
}


print("Combining the bootstrap outs:")


path_out_boot <- "results_bootstrap"
setwd(path_out_boot)

for (group_oi in names(type_groups)){
  
  print(group_oi)
  
  files_boot <- dir(pattern = glob2rx(type_groups[group_oi]))
  
  df_oi <- read.table(files_boot[1],header=TRUE)
  n_rows <- dim(df_oi)[1]
  
  df_boot <- data.frame(matrix(NA, nrow = length(files_boot)*n_rows, ncol = dim(df_oi)[2]))
  colnames(df_boot) <- colnames(df_oi)
  
  counter <- 1
  problem_files <- c()
  for (file_oi in files_boot){
    print(file_oi)
    
    df_oi <- read.table(file_oi,header=TRUE)
    if (dim(df_oi)[1] != n_rows){
      print(sprintf("Something weird with: %s",file_oi))
      problem_files <- c(problem_files,file_oi)
    } else{
      idx <- seq(n_rows)+(counter-1)*n_rows
      df_boot[idx,] <- df_oi
    }
    counter <- counter+1
  }
  
  saveRDS(df_boot,sprintf("table_bootstrap_minP_noP_all_CREs_%s_all_reps_v3_20220823.RDS",group_oi))
  
  
}


