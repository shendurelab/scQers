

file_oBC_mBC <- "/Users/jbl/Documents/UW/manuscripts/single_cell_expression_reporter/GEO_processed_files/p025_recloned_complex_mBC_oBC_subassembly.txt"
df_oBC_mBC_p25 <- read.table(file_oBC_mBC,sep="\t",header=TRUE)
df_oBC_mBC_p25_valid <- df_oBC_mBC_p25 %>% dplyr::filter(bool_valid) %>% select(mBC,oBC)



setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE")

# sample_name <- "chimeric_CRE_A_S6"
# sample_name <- "chimeric_CRE_B_S5"
# sample_name <- "chimeric_CRE_C_S4"
# sample_name <- "chimeric_CRE_D_S3"


file_in <- sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq037_chimeric_and_mutated_CRE_subassembly_20230218_NS2000/chimericCRE/pileup_run2/intermediate/%s_read2_pileup_k2_20230301.txt.gz",
                   sample_name)
df_oBC_CRE1 <- read.table(file_in,header=TRUE)
  
  
  
count_thresh <- 5

dev.set(4)
ggplot(df_oBC_CRE1) + stat_bin(aes(x=n_counts),geom="step")+
  geom_vline(xintercept = count_thresh+1)+
  scale_x_log10()+scale_y_log10()

df_oBC_CRE1_hi <- df_oBC_CRE1 %>% filter(n_counts>count_thresh) %>% rename(BC="RC_oBC")

# dim(df_oBC_CRE1_hi)
# 

# head(df_oBC_CRE1_hi)

df_oBC_CRE1_hi$oBC <-  DNAStringSet(df_oBC_CRE1_hi$RC_oBC) %>% reverseComplement() %>% as.character()


mean(df_oBC_CRE1_hi$oBC %in% df_oBC_mBC_p25_valid$oBC)
mean(df_oBC_CRE1_hi$oBC %in% df_oBC_mBC_p25$oBC)

df_oBC_CRE1_hi2 <- df_oBC_CRE1_hi %>% transform(non_unique_oBC_CRE=(duplicated(oBC) | duplicated(oBC, fromLast = TRUE)))
mean(df_oBC_CRE1_hi2$non_unique_oBC_CRE)

df_oBC_CRE1_hi2 %>% filter(non_unique_oBC_CRE) %>% arrange(oBC)


df_oBC_CRE1_hi3 <- df_oBC_CRE1_hi2 %>% left_join(df_oBC_mBC_p25 %>% select(oBC,mBC,valid_unique_oBC_mBC_pair=bool_valid) )
  
length(df_oBC_CRE1_hi$RC_oBC %>% unique)/length(df_oBC_CRE1_hi$RC_oBC)

write.table(df_oBC_CRE1_hi3,sprintf("oBC_CRE1_hi_counts_w_paired_mBC_%s_20230820.txt",sample_name),
            row.names=FALSE,quote=FALSE,sep="\t")
  







# # # # # # # # # # # # # # # # # # # # # # # # 
# hamming distance to nanopore BCs?
# # # # # # # # # # # # # # # # # # # # # # # # 

# sample_nanop <- "BC1_chimericCRE"
# file_in_illumina <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE/oBC_CRE1_hi_counts_w_paired_mBC_chimeric_CRE_A_S6_20230820.txt"

# sample_nanop <- "BC2_chimericCRE"
# file_in_illumina <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE/oBC_CRE1_hi_counts_w_paired_mBC_chimeric_CRE_C_S4_20230820.txt"


# sample_nanop <- "BC3_chimericCRE"
# file_in_illumina <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE/oBC_CRE1_hi_counts_w_paired_mBC_chimeric_CRE_B_S5_20230820.txt"
# file_nanop_in <- "/Users/jbl/Documents/UW/Data/nanopore/nano01_chimeric_CREs_20230523/nanop001_chimCRE_BC3_mBC_CRE_pileup_20230702.txt"

sample_nanop <- "BC4_chimericCRE"
file_in_illumina <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE/oBC_CRE1_hi_counts_w_paired_mBC_chimeric_CRE_D_S3_20230820.txt"
file_nanop_in <- "/Users/jbl/Documents/UW/Data/nanopore/nano01_chimeric_CREs_20230523/nanop001_chimCRE_BC4_mBC_CRE_pileup_20230702.txt"


df_oBC_CRE_mBC_illumina <- read.table(file_in_illumina,header=TRUE)
df_nanop_pile <- read.table(file_nanop_in,header=TRUE)



dim(df_nanop_pile %>% filter(mBC_counts>=2))

mean(df_nanop_pile$mBC %in% df_oBC_CRE_mBC_illumina$mBC)

# min_lv_dist_BCs <- vector(mode="numeric",length=dim(df_nanop_pile)[1])

mBC_illumina <- df_oBC_CRE_mBC_illumina %>% filter(!is.na(mBC)) %>% pull(mBC) %>% unique

df_nanop_pile2 <- df_nanop_pile
for (ind1 in seq(1,dim(df_nanop_pile)[1])){
  
  lv_dist_illumina_mBC <- stringdist(df_nanop_pile2$mBC[ind1],mBC_illumina,method="lv")
  
  df_nanop_pile2$min_lv_dist_BC[ind1] <- min(lv_dist_illumina_mBC)
  df_nanop_pile2$n_min_dist[ind1] <- sum(lv_dist_illumina_mBC==min(lv_dist_illumina_mBC))
  
  if (sum(lv_dist_illumina_mBC==min(lv_dist_illumina_mBC))==1){
    df_nanop_pile2$error_corr_mBC[ind1]<- mBC_illumina[min(lv_dist_illumina_mBC)==lv_dist_illumina_mBC]
  } else {
    df_nanop_pile2$error_corr_mBC[ind1] <- NA
  }
  
  if ((ind1 %% 500)==0){
    print(ind1)
  }
}

write.table(df_nanop_pile2,sprintf("df_nanop_pile2_%s_20230820.txt",sample_nanop),
            row.names=FALSE,quote=FALSE,sep="\t")



library(universalmotif)


df_shuffled_seq <- df_nanop_pile

for (ind1 in seq(1,1000)){
  
  lv_dist_illumina_mBC <- stringdist(shuffle_sequences(DNAString(df_nanop_pile$mBC[ind1])) %>% as.character,mBC_illumina,method="lv")
  
  df_shuffled_seq$min_lv_dist_BC[ind1] <- min(lv_dist_illumina_mBC)
  df_shuffled_seq$n_min_dist[ind1] <- sum(lv_dist_illumina_mBC==min(lv_dist_illumina_mBC))
  
  if (sum(lv_dist_illumina_mBC==min(lv_dist_illumina_mBC))==1){
    df_shuffled_seq$error_corr_mBC[ind1]<- mBC_illumina[min(lv_dist_illumina_mBC)==lv_dist_illumina_mBC]
  } else {
    df_shuffled_seq$error_corr_mBC[ind1] <- NA
  }
  
  
  if ((ind1 %% 500)==0){
    print(ind1)
  }
  
}




# df_nanop_pile2 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq047_reseq_scv2_lenti_scQer_Tony_MPRA_Hanna_sci_20230728/misc_subassembly_chimeric_CRE/df_nanop_pile2_BC2_chimericCRE_20230820.txt",header=TRUE)


dev.set(6)
# ggplot(df_nanop_pile2) + stat_bin2d(aes(x=mBC_counts,y=min_lv_dist_BC))+
#   scale_x_log10()

ggplot(df_shuffled_seq[1:1000,]) + geom_step(aes(x=min_lv_dist_BC),stat="ecdf")

df_nanop_pile2$count_cat <- "none"
df_nanop_pile2$count_cat[df_nanop_pile2$mBC_counts<2] <- "low_1"
df_nanop_pile2$count_cat[df_nanop_pile2$mBC_counts>=2 & df_nanop_pile2$mBC_counts<4] <- "mid_2-3"
df_nanop_pile2$count_cat[df_nanop_pile2$mBC_counts>=4] <- "hi_4+"

# ggplot(df_nanop_pile2) + stat_bin(aes(x=min_lv_dist_BC),geom="step")
ggplot(df_nanop_pile2) + geom_step(aes(x=min_lv_dist_BC,color=count_cat),stat="ecdf")





min_lv_dist_thresh <- 2
read_thresh <- 1

df_nanop_pile3 <- df_nanop_pile2 %>% filter(min_lv_dist_BC<=min_lv_dist_thresh & !is.na(error_corr_mBC) & mBC_counts>read_thresh) %>%
  left_join(df_oBC_CRE_mBC_illumina %>% filter(!non_unique_oBC_CRE) %>%
              select(error_corr_mBC=mBC,
                     CRE_id_ill=CRE_id,
                     CRE_ori_ill=read_forward_in_ref,
                     valid_unique_oBC_mBC_pair))

df_nanop_pile4 <- df_nanop_pile3 %>% transform(CRE_ori_ill2=ifelse(CRE_ori_ill=="True",'-','+')) %>% select(-CRE_ori_ill)

# str_content <- df_nanop_pile4$CRE_pair_id %>% str_split("::")
# CRE_L_w_ori <- str_content %>% lapply("[[",1) %>% unlist()
# CRE_L_id <- CRE_L_w_ori %>% substr(1,str_length(CRE_L_w_ori)-1)
# CRE_L_ori <- CRE_L_w_ori %>% substr(str_length(CRE_L_w_ori),str_length(CRE_L_w_ori))
# CRE_R_w_ori <- str_content %>% lapply("[[",2) %>% unlist()
# CRE_R_id <- CRE_R_w_ori %>% substr(1,str_length(CRE_R_w_ori)-1)
# CRE_R_ori <- CRE_R_w_ori %>% substr(str_length(CRE_R_w_ori),str_length(CRE_R_w_ori))
# df_nanop_pile4$CRE_L_id <- CRE_L_id
# df_nanop_pile4$CRE_L_ori <- CRE_L_ori
# df_nanop_pile4$CRE_R_id <- CRE_R_id
# df_nanop_pile4$CRE_R_ori <- CRE_R_ori


mean(df_nanop_pile4$CRE_L_id==df_nanop_pile4$CRE_id_ill,na.rm=TRUE)
mean(df_nanop_pile4$CRE_L_ori==df_nanop_pile4$CRE_ori_ill2,na.rm=TRUE)

df_nanop_pile_err_corr <- df_nanop_pile4 %>% filter(CRE_L_id==CRE_id_ill & CRE_L_ori==CRE_ori_ill2)

table(df_nanop_pile_err_corr$count_cat)

df_nanop_pile_err_corr2 <- df_nanop_pile_err_corr %>% select(-c(mBC,min_lv_dist_BC,mBC_counts,n_min_dist,count_cat)) %>% unique() %>% 
  arrange(error_corr_mBC) %>% transform(non_unique_err_corr_mBC=(duplicated(error_corr_mBC) | duplicated(error_corr_mBC, fromLast = TRUE))) %>% 
  filter(valid_unique_oBC_mBC_pair & !non_unique_err_corr_mBC) %>% select(-c(CRE_id_ill,valid_unique_oBC_mBC_pair,CRE_ori_ill2,non_unique_err_corr_mBC))

dim(df_nanop_pile_err_corr2)

# df_nanop_pile_err_corr2_lv2 <- df_nanop_pile_err_corr2
# df_nanop_pile_err_corr2_lv1 <- df_nanop_pile_err_corr2
# 
# head(df_nanop_pile_err_corr2_lv2)
# 
# mean(df_nanop_pile_err_corr2_lv2$error_corr_mBC %in% df_nanop_pile_err_corr2_lv1$error_corr_mBC)
# mean(df_nanop_pile_err_corr2_lv1$error_corr_mBC %in% df_nanop_pile_err_corr2_lv2$error_corr_mBC)



write.table(df_nanop_pile_err_corr2,
            sprintf("error_corrected_nanopore_subassembly_readThresh%d_lvDist%d_%s_20230820.txt",read_thresh,min_lv_dist_thresh,sample_nanop),quote=FALSE,row.names=FALSE,sep="\t")

# df_suba <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq046_scQer_scv2_mEB_10x_20230713/final_subassembly_scv2_no_collision_unique_matched_mBC_oBC_20230714.txt",
#                       header=TRUE)
# table(df_suba$lib_id)

df_nanop_pile_err_corr2[c(19,20),]
