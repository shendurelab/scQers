#library(dplyr)
#library(patchwork)
#library(sctransform)
#library(ggplot2)
library(tidyverse)
#library(scales)
library(tictoc)
library(stringdist)
library(igraph)



args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]
file_name_corrected <- args[2]
variable_name <- args[3]

# file_name <- "/Users/jbl/Documents/UW/Data/sequencing_runs/20210828_10x_pigglyFlex_MOI_validation/valid_cBC_get_bc_files/L1_get_bc_v2_poly_dT_valid_cBC_20211006.txt"
# file_name_corrected <- "/Users/jbl/Documents/UW/Data/sequencing_runs/20210828_10x_pigglyFlex_MOI_validation/valid_cBC_get_bc_files/L1_get_bc_v2_poly_dT_valid_cBC_corrected_UMIs_20211006.txt"
# variable_name <- "mBC"

df_get_bc_v2 <- read.table(file_name,header=TRUE, sep="\t")

print(head(df_get_bc_v2))

# function with Hamming distance calculation and connected component calculation
conn_components_UMI <-function(UMI_str,dist_thresh){
  
  # parse the string
  df_UMI <- data.frame(umi_reads=str_split(UMI_str,",")[[1]])
  df_UMI$umi <- sapply(str_split(df_UMI$umi_reads,"_"),"[")[1,]
  df_UMI$count <- strtoi(sapply(str_split(df_UMI$umi_reads,"_"),"[")[2,])
  umis <- df_UMI$umi
  
  # matrix of hamming distance between UMIs (VECTORIZED)
  Ham_dist_UMIs <- matrix(data=0,nrow=length(umis),ncol=length(umis))
  for (ind1 in seq(length(umis))){
    Ham_dist_UMIs[ind1,] <- stringdist(umis[ind1],umis,method="hamming")
  }
  
  # rownames(Ham_dist_UMIs) <- umis
  # colnames(Ham_dist_UMIs) <- umis
  # for (ind1 in seq(length(umis))){
  #   for (ind2 in seq(length(umis))){
  #     if (ind2>ind1){
  #       umi1 <- umis[ind1]
  #       umi2 <- umis[ind2]
  #       Ham_dist_UMIs[umi1,umi2] <- stringdist(umi1,umi2,method="hamming")
  #     }
  #   }
  # }
  # Ham_dist_UMIs[lower.tri(Ham_dist_UMIs)] = t(Ham_dist_UMIs)[lower.tri(Ham_dist_UMIs)]
  
  # connected components
  conn_mat <- Ham_dist_UMIs<=dist_thresh
  g  <- graph.adjacency(conn_mat)
  clu <- components(g)
  return(clu$no)
}




dist_thresh <- 1
n_cBC_mBC <- dim(df_get_bc_v2)[1]

corrected_UMI <- vector(mode="numeric",length=n_cBC_mBC)

for (ind in seq(n_cBC_mBC)){
  
  if (df_get_bc_v2$n_UMI_filtered[ind]>1){
    filtered_UMIs_oi <- df_get_bc_v2$list_reads_per_UMIs_filtered[ind]
    corrected_UMI[ind] <- conn_components_UMI(filtered_UMIs_oi,dist_thresh)
  } else {
    corrected_UMI[ind] <- 1
  }
  
  if ((ind %% ceiling(0.001*n_cBC_mBC))==0){
    print(ind/n_cBC_mBC)
  }
}

# adding corrected UMIs to df and writing output to text file
df_get_bc_v2_2 <- df_get_bc_v2 %>% transform(filtered_corrected_UMIs=corrected_UMI)
df_out <- df_get_bc_v2_2 %>% select(c("cBC",variable_name,"n_reads_filtered","n_UMI_filtered","filtered_corrected_UMIs"))
write.table(df_out,file_name_corrected, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")






