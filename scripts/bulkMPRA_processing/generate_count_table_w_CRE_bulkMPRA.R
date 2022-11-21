require(dplyr); require(data.table); require(ggplot2);
require(cowplot); theme_set(theme_cowplot()); require(stringr)
library(scales)
library(Biostrings)
library(patchwork)
library(stringdist)
library(scales)
library(tidyverse)
library(tictoc)

library(stringr)

# start processing the list of barcodes

dir_out <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq032_bottlenecked_mEB_bulk_MPRA_rep2B_20220613/rep2B_MPRA_outs_20220614/"
setwd(dir_out)

dir_in <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq032_bottlenecked_mEB_bulk_MPRA_rep2B_20220613/pre_processing_mEB_MPRA_20220614/BC_quant/"
files <- dir(dir_in)

sample_name <- str_replace(str_replace(files,"bulk_MPRA_mEBs_",""),"_BC_quant_20220614.txt.gz","") 
names(sample_name) <- files

splt_sample_name <- sample_name %>% str_split("_") 

sample_type <- splt_sample_name %>% lapply("[[",3) %>% unlist
rep_id <- splt_sample_name %>% lapply("[[",1) %>% unlist %>% substr(4,6)
time_point <- splt_sample_name %>% lapply("[[",2) %>% unlist

prep_batch <- rep("batch3",length(sample_type))

df_sample_data <- data.frame(sample_name,sample_type,rep_id,time_point,prep_batch)
rownames(df_sample_data) <- sample_name



# list of final barcodes
file_suba <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq023_devCRE_subassemb_sci_mBC_oBC_MPRA_20220331/oBC_devCRE_subassembly_final_outs_w_multimappers_20220818/final_list_v2_subassembly_oBC_CRE_mBC_triplets_20220818.txt"
df_mBC_pool2 <- read.table(file_suba, header=TRUE) %>% select(mBC,CRE_id,CRE_class)
date_str <- "20220826"


df_MPRA <- data.frame()


for (file_oi in files){
  
  print(file_oi)
  
  
  df_BC_oi <- read.table(paste0(dir_in,file_oi),header=TRUE)
  
  
  ## QC plots and assessment
  thresh_reads <- 10
  df_BC_oi2 <- df_BC_oi %>% filter(n_reads_per_mBC>thresh_reads)
  
  print(sample_name[file_oi] %>% as.character())

  print("Number of mBC >10 reads:")
  print(dim(df_BC_oi2)[1])
  
  print("proportion mBC (>10 reads) devCRE:")
  print(mean(df_BC_oi2$mBC %in% (df_mBC_pool2 %>% filter(CRE_class=="devCRE") %>% pull(mBC) )))
  
  print("proportion mBC (>10 reads) promoter:")
  print(mean(df_BC_oi2$mBC %in% (df_mBC_pool2 %>% filter(CRE_class=="promoters") %>% pull(mBC) )))

  
  df_BC_oi3 <- df_BC_oi2 %>% left_join(df_mBC_pool2)
  
  
 
  # plt_cdf <- ggplot(df_BC_oi)+geom_step(aes(x=n_reads_per_mBC,y=1-..y..),stat="ecdf")+scale_x_log10()+scale_y_log10()+
  #   labs(title=sample_name[file_oi])
  # 
  # 
  # fig_name <- sprintf("%sfig_cdf_mBC_counts_%s_20220614.pdf",dir_out,sample_name[file_oi])
  # pdf(fig_name,width=7,height=5)
  # print(plt_cdf)
  # dev.off()
  # 
  # plt_reads_per_BC <- ggplot(df_BC_oi3)+stat_bin(aes(x=n_reads_per_mBC,color=CRE_class),
  #                                position="identity",geom="step",bins=50)+scale_x_log10()+scale_y_log10()+
  #   labs(title=sample_name[file_oi])
  # 
  # fig_name <- sprintf("%sfig_dist_mBC_by_class_%s_20220614.pdf",dir_out,sample_name[file_oi])
  # pdf(fig_name,width=7,height=5)
  # print(plt_reads_per_BC)
  # dev.off()
  
 
  names(df_BC_oi) <- c("mBC","umi","reads")
  df_MPRA_oi<- df_mBC_pool2 %>% left_join(df_BC_oi)
  
  df_MPRA_oi2 <- df_MPRA_oi %>% transform(norm_umi=umi/sum(umi,na.rm=TRUE),norm_reads=reads/sum(reads,na.rm=TRUE))
  
  df_sample_metadata_oi <- df_sample_data[sample_name[file_oi],] %>% dplyr::slice(rep(1:n(), each = dim(df_MPRA_oi)[1]))
  
  df_MPRA_oi3 <- df_MPRA_oi2 %>% cbind(df_sample_metadata_oi)

  df_MPRA <- rbind(df_MPRA,df_MPRA_oi3)
  
  
}


# try to pivot wider for the RNA and DNA

df_MPRA2 <- df_MPRA %>% select(-sample_name)

df_MPRA_wide <- df_MPRA2 %>% 
  pivot_wider(names_from=sample_type,values_from=c(umi,reads,norm_umi,norm_reads))

write.table(df_MPRA_wide,sprintf("bulk_MPRA_mEB_rep2B_batch3_seq032_v2_%s.txt",date_str),
            sep="\t",quote=FALSE,row.names=FALSE)

