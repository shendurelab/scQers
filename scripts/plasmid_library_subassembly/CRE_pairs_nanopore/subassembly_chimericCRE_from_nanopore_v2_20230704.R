
library(tidyverse)
library(stringr)
library(ShortRead)
library(seqinr)
library(tictoc)
library(stringdist)
library(igraph)
library(data.table)



# # # # # # # # # # # # # # # # # # # # # # # #
# input files and parameters 
# # # # # # # # # # # # # # # # # # # # # # # #

signpost_score_thresh <- 50 # for region of 30, max score is 30*2 = 60
n_signpost_mapped_thresh <- 4 # we have 5 signposts

size_lo_CRE <- 500
size_hi_CRE <- 1500
size_lo_BC <- 42
size_hi_BC <- 52

aln_CRE_thresh <- 750



setwd("/Users/jbl/Documents/UW/Data/nanopore/nano01_chimeric_CREs_20230523/")
signpost_fasta <- "signposts_chimericCREs_20230627.fasta"
CRE_fasta <- "chimeric_CRE_components_20230701.fasta"
date_str <- "20230704"


input_fastq <- "nanop001_chimCRE_BC1_20230627.fastq.gz"
output_file_id <- "nanop001_BC1_chimericCRE" 

input_fastq <- "nanop001_chimCRE_BC2_20230627.fastq.gz"
output_file_id <- "nanop001_BC2_chimericCRE" 



if (TRUE){

# # # # # # # # # # # # # # # # # # # # # # # #
# upstream signpost aln w/ SW
# # # # # # # # # # # # # # # # # # # # # # # #
file_aln_signposts <- sprintf("%s_signposts_aln_fmt2_%s.txt",output_file_id,date_str)

bsh_cmd <- sprintf("./sw_aln_reformat_20230704.sh %s %s %s",input_fastq,signpost_fasta,file_aln_signposts)
system(bsh_cmd,wait=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# load and initial processing of signpost alignment file
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_aln <- read.table(file_aln_signposts,sep="\t")
df_aln2 <- df_aln
for (col_id in seq(1,dim(df_aln)[2])){
  
  col_content <- str_split(df_aln[,col_id],": ")
  
  col_name <- col_content[[1]][1]
  col_var_content <- col_content %>% lapply("[[",2) %>% unlist
  
  num_col_var_content <- as.numeric(col_var_content)
  
  if (sum(is.na(num_col_var_content))==0){
    df_aln2[,col_id] <- num_col_var_content
  } else {
    df_aln2[,col_id] <- col_var_content
  }
  colnames(df_aln2)[col_id] <- col_name
  
}

df_aln2 %>% pull(target_name) %>% unique %>% length

# optional plots
# ggplot(df_aln2) + 
#   stat_bin2d(aes(y=optimal_alignment_score,x=suboptimal_alignment_score,fill=log(..count..)))+
#   facet_wrap(~query_name,nrow=1)+
#   coord_fixed()
# 
# 
# ggplot(df_aln2) + 
#   stat_bin(aes(x=optimal_alignment_score),bins=60)+
#   facet_wrap(~query_name,nrow=1)


# identifying reads with high alignment score for multiple signposts
df_aln2_hi <- df_aln2 %>% filter(optimal_alignment_score>=signpost_score_thresh) %>% arrange(target_name)
df_reads <- df_aln2_hi %>% group_by(target_name) %>% summarize(n_signpost=length(query_name))
df_aln2_hi2 <- df_aln2 %>% filter(target_name %in% (df_reads %>% filter(n_signpost>=n_signpost_mapped_thresh) %>% pull(target_name))) %>% arrange(target_name)
df_aln2_hi2 %>% pull(target_name) %>% unique %>% length


# stratifying by various signposts
df_aln2_hi2_BC_L <- df_aln2_hi2 %>% filter(query_name=="BC_L") %>% select(target_name,pos_BC_L=target_end)
df_aln2_hi2_BC_R <- df_aln2_hi2 %>% filter(query_name=="BC_R") %>% select(target_name,pos_BC_R=target_end)
df_aln2_hi2_CRE_L <- df_aln2_hi2 %>% filter(query_name=="CRE_L") %>% select(target_name,pos_CRE_L=target_end)
df_aln2_hi2_CRE_R <- df_aln2_hi2 %>% filter(query_name=="CRE_R") %>% select(target_name,pos_CRE_R=target_end)
df_aln2_hi2_CRE_J <- df_aln2_hi2 %>% filter(query_name=="CRE_J") %>% select(target_name,pos_CRE_J=target_end, strand)

df_aln2_hi_pos <- df_aln2_hi2_BC_L %>% 
  left_join(df_aln2_hi2_BC_R) %>%
  left_join(df_aln2_hi2_CRE_L) %>%
  left_join(df_aln2_hi2_CRE_R) %>%
  left_join(df_aln2_hi2_CRE_J)
df_aln2_hi_pos2 <- df_aln2_hi_pos %>% transform(len_BC=abs(pos_BC_R-pos_BC_L),
                                                len_CRE_L=abs(pos_CRE_L-pos_CRE_J),
                                                len_CRE_R=abs(pos_CRE_J-pos_CRE_R))

# ggplot(df_aln2_hi_pos2) + stat_bin(aes(x=len_BC))+xlim(c(35,55))
# 
# ggplot(df_aln2_hi_pos2) + 
#   stat_bin(aes(x=len_CRE_L),color="red",geom="step")+xlim(c(0,2000))+
#   stat_bin(aes(x=len_CRE_R),color="blue",geom="step")+xlim(c(0,2000))


# selecting reads with roughly correctly sized signposts
df_usable_reads <- df_aln2_hi_pos2 %>% 
  filter(len_BC>size_lo_BC & len_BC<size_hi_BC) %>%
  filter(len_CRE_R>size_lo_CRE & len_CRE_R<size_hi_CRE) %>%
  filter(len_CRE_L>size_lo_CRE & len_CRE_L<size_hi_CRE)



# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# extract BC and CRE sequences from fastq raw reads
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_fastq <- readFastq(input_fastq)

nanop_read_l <- str_length(df_fastq@sread)
# ggplot()+stat_bin(aes(x=nanop_read_l),geom="step")+scale_x_log10()+annotation_logticks()
# sum(nanop_read_l>5000 & nanop_read_l<10000)

df_usable_reads_w_seq <- df_usable_reads
df_usable_reads_w_seq$seq_CRE_L <- NA
df_usable_reads_w_seq$seq_CRE_R <- NA
df_usable_reads_w_seq$seq_BC <- NA


df_fastq_usable <- df_fastq[df_fastq@id %in% df_usable_reads$target_name]
valid_read_id <- df_fastq_usable@id %>% as.character

# parsing fastq sequences with identifiable signposts, extract relevant sequences 
counter <- 0
l_signpost <- 30
for (read_id in valid_read_id){
  
  
  if (counter%%1000==0){
    print(sprintf("%s, %d/%d=%.3f",read_id,counter,length(valid_read_id),counter/length(valid_read_id)))
  }
  
  df_usable_reads_oi <- df_usable_reads %>% filter(target_name==read_id)
  interval_CRE_L <- c(df_usable_reads_oi$pos_CRE_L,df_usable_reads_oi$pos_CRE_J) %>% sort()
  interval_CRE_R <- c(df_usable_reads_oi$pos_CRE_J,df_usable_reads_oi$pos_CRE_R) %>% sort()
  interval_BC <- c(df_usable_reads_oi$pos_BC_L,df_usable_reads_oi$pos_BC_R) %>% sort()
  
  
  # get sequence
  df_fastq_oi <- df_fastq_usable[df_fastq_usable@id==read_id]
  read_seq <- df_fastq_oi@sread
  
  # get reverse complement if reverse orientation, and correct for offset
  if (df_usable_reads_oi$strand=="-"){
    seq_CRE_L_oi <- reverseComplement(subseq(read_seq, start=max(interval_CRE_L[1]-l_signpost+1,1), end=interval_CRE_L[2]-l_signpost+1))
    seq_CRE_R_oi <- reverseComplement(subseq(read_seq, start=max(interval_CRE_R[1]-l_signpost+1,1), end=interval_CRE_R[2]-l_signpost+1))
    seq_BC_oi <- reverseComplement(subseq(read_seq, start=max(interval_BC[1]-l_signpost+1,1), end=interval_BC[2]-l_signpost+1))
  } else {
    seq_CRE_L_oi <- subseq(read_seq, start=interval_CRE_L[1], end=interval_CRE_L[2])
    seq_CRE_R_oi <- subseq(read_seq, start=interval_CRE_R[1], end=interval_CRE_R[2])
    seq_BC_oi <- subseq(read_seq, start=interval_BC[1], end=interval_BC[2])
  }
    
  
  # assign
  idx <- which(df_usable_reads_w_seq$target_name==read_id)
  df_usable_reads_w_seq[idx,]$seq_CRE_L <- seq_CRE_L_oi
  df_usable_reads_w_seq[idx,]$seq_CRE_R <- seq_CRE_R_oi
  df_usable_reads_w_seq[idx,]$seq_BC <- seq_BC_oi
  
  counter <- counter+1
  
}


# print intermediate files: 
write.table(df_usable_reads_w_seq,
            sprintf("%s_usable_reads_w_seq_%s.txt",output_file_id,date_str),
            sep="\t",quote=FALSE,row.names=FALSE)

write.table(df_usable_reads_w_seq %>% select(target_name,strand,seq_BC),
            sprintf("%s_usable_reads_w_seq_BC_only_%s.txt",output_file_id,date_str),
            sep="\t",quote=FALSE,row.names=FALSE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# create fasta file of CRE seq and align to expected CREs
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

file_CRE_L_fasta <- sprintf("%s_CRE_L_%s.fasta",output_file_id,date_str)
file_CRE_R_fasta <- sprintf("%s_CRE_R_%s.fasta",output_file_id,date_str)

write.fasta(as.list(df_usable_reads_w_seq$seq_CRE_L),
            names=df_usable_reads_w_seq$target_name,
            file.out=file_CRE_L_fasta)

write.fasta(as.list(df_usable_reads_w_seq$seq_CRE_R),
            names=df_usable_reads_w_seq$target_name,
            file.out=file_CRE_R_fasta)


file_aln_CRE_L <- sprintf("%s_CRE_L_aln_fmt_%s.txt",output_file_id,date_str)
file_aln_CRE_R <- sprintf("%s_CRE_R_aln_fmt_%s.txt",output_file_id,date_str)

bsh_cmd_CRE_L <- sprintf("./sw_aln_reformat_20230704.sh %s %s %s",file_CRE_L_fasta,CRE_fasta,file_aln_CRE_L)
system(bsh_cmd_CRE_L,wait=TRUE)

bsh_cmd_CRE_R <- sprintf("./sw_aln_reformat_20230704.sh %s %s %s",file_CRE_R_fasta,CRE_fasta,file_aln_CRE_R)
system(bsh_cmd_CRE_R,wait=TRUE)




# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# final CRE assignment based on SW alignment
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # #  CRE_R
df_CRE_R_aln <- read.table(file_aln_CRE_R,header=FALSE,sep="\t")

# reformat
df_CRE_R_aln2 <- df_CRE_R_aln
for (col_id in seq(1,dim(df_CRE_R_aln)[2])){
  
  col_content <- str_split(df_CRE_R_aln[,col_id],": ")
  col_name <- col_content[[1]][1]
  col_var_content <- col_content %>% lapply("[[",2) %>% unlist
  num_col_var_content <- as.numeric(col_var_content)
  
  if (sum(is.na(num_col_var_content))==0){
    df_CRE_R_aln2[,col_id] <- num_col_var_content
  } else {
    df_CRE_R_aln2[,col_id] <- col_var_content
  }
  colnames(df_CRE_R_aln2)[col_id] <- col_name
}

df_CRE_R_aln3 <- df_CRE_R_aln2 %>% arrange(target_name)
ggplot(df_CRE_R_aln3) + stat_bin(aes(x=optimal_alignment_score),geom="step")+scale_y_log10()

df_max_aln_CRE_R <- df_CRE_R_aln3 %>% group_by(target_name) %>% slice_max(optimal_alignment_score)
df_max_aln_CRE_R2 <- df_max_aln_CRE_R
df_max_aln_CRE_R2[df_max_aln_CRE_R2$optimal_alignment_score<aln_CRE_thresh,"query_name"] <- NA

ggplot(df_max_aln_CRE_R2) + stat_bin(aes(x=optimal_alignment_score,color=query_name),geom="step",position="identity")+scale_y_log10()

df_max_aln_CRE_R <- df_max_aln_CRE_R2 %>% select(target_name, CRE_R_id=query_name, CRE_R_ori=strand)

table(df_max_aln_CRE_R$CRE_R_id, df_max_aln_CRE_R$CRE_R_ori)


# # # # # # # # # #  CRE_L
df_CRE_L_aln <- read.table(file_aln_CRE_L,header=FALSE,sep="\t")

# reformat
df_CRE_L_aln2 <- df_CRE_L_aln
for (col_id in seq(1,dim(df_CRE_L_aln)[2])){
  
  col_content <- str_split(df_CRE_L_aln[,col_id],": ")
  col_name <- col_content[[1]][1]
  col_var_content <- col_content %>% lapply("[[",2) %>% unlist
  num_col_var_content <- as.numeric(col_var_content)
  
  if (sum(is.na(num_col_var_content))==0){
    df_CRE_L_aln2[,col_id] <- num_col_var_content
  } else {
    df_CRE_L_aln2[,col_id] <- col_var_content
  }
  colnames(df_CRE_L_aln2)[col_id] <- col_name
}

df_CRE_L_aln3 <- df_CRE_L_aln2 %>% arrange(target_name)
ggplot(df_CRE_L_aln3) + stat_bin(aes(x=optimal_alignment_score),geom="step")#+scale_y_log10()

df_max_aln_CRE_L <- df_CRE_L_aln3 %>% group_by(target_name) %>% slice_max(optimal_alignment_score)
df_max_aln_CRE_L2 <- df_max_aln_CRE_L
df_max_aln_CRE_L2[df_max_aln_CRE_L2$optimal_alignment_score<aln_CRE_thresh,"query_name"] <- NA

ggplot(df_max_aln_CRE_L2) + stat_bin(aes(x=optimal_alignment_score,color=query_name),geom="step",position="identity")+scale_y_log10()

df_max_aln_CRE_L <- df_max_aln_CRE_L2 %>% select(target_name, CRE_L_id=query_name, CRE_L_ori=strand)

table(df_max_aln_CRE_L$CRE_L_id, df_max_aln_CRE_L$CRE_L_ori)



# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# join with original table
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# df_usable_reads_w_seq2 <- read.table("nanop001_chimCRE_BC4_usable_reads_w_seq_20230702.txt",
#             sep="\t",header=TRUE)

df_usable_reads_w_seq_mini <- df_usable_reads_w_seq %>% select(target_name,seq_BC)

# extract the actual mBC sequence & join CRE info
df_mBC_w_CRE <- df_usable_reads_w_seq_mini %>% transform(mBC=substr(seq_BC,2,16)) %>% select(-seq_BC) %>% 
  left_join(df_max_aln_CRE_L) %>% left_join(df_max_aln_CRE_R)

df_mBC_w_CRE2 <- df_mBC_w_CRE %>% transform(CRE_pair_id=paste0(CRE_L_id,CRE_L_ori,"::",CRE_R_id,CRE_R_ori))




df_mBC_w_CRE2_pile <- df_mBC_w_CRE2 %>% group_by(CRE_L_id,CRE_L_ori,CRE_R_id,CRE_R_ori,CRE_pair_id,mBC) %>% summarize(mBC_counts=length(mBC))

file_mBC_CRE_pileup <- sprintf("%s_mBC_CRE_pileup_%s.txt",output_file_id,date_str)
write.table(df_mBC_w_CRE2_pile %>% arrange(mBC),file_mBC_CRE_pileup,
            sep="\t",quote=FALSE,row.names=FALSE)
}

# df_mBC_w_CRE2_pile2 <- df_mBC_w_CRE2_pile %>% filter(CRE_L_id != CRE_R_id)

# ggplot() +
#   stat_bin(aes(x=df_mBC_w_CRE2_pile$mBC_counts),geom="step",color="black")+scale_x_log10()+scale_y_log10()
  # geom_vline(aes(xintercept=mBC_max_read_thresh),linetype="dotted")



# # build levenshtein communities
# 
# CRE_pairs <- df_mBC_w_CRE2_pile2 %>% pull(CRE_pair_id) %>% unique
# pair_id <- 1
# 
# df_list <- list()
# for (CRE_pair_oi in CRE_pairs){
#   
#   print(CRE_pair_oi)
#   df_mBC_w_CRE2_pile2_oi <- df_mBC_w_CRE2_pile2 %>% filter(CRE_pair_id==CRE_pair_oi)
#   
#   mBCs <- df_mBC_w_CRE2_pile2_oi$mBC
#   dist_BCs <- matrix(data=0,nrow=length(mBCs),ncol=length(mBCs))
#   
#   for (ind1 in seq(length(mBCs))){
#     
#     # hamm_BCs[ind1,] <- stringdist(paired_BCs[ind1],paired_BCs,method="hamming")
#     dist_BCs[ind1,] <- stringdist(mBCs[ind1],mBCs,method="lv")
#     
#     if ((ind1 %% 500)==0){
#       print(ind1)
#     }
#   }
#   
#   dist_thresh <- 2
#   conn_mat <- dist_BCs<=dist_thresh
#   g  <- graph.adjacency(conn_mat)
#   clu <- components(g)
#   
#   # adding membership
#   df_mBC_w_CRE2_pile2_oi2 <- df_mBC_w_CRE2_pile2_oi %>% transform(conn_clu_id=sprintf("%d_%d",pair_id,clu$membership)) %>% arrange(conn_clu_id)
#   
#   df_list[[pair_id]] <- df_mBC_w_CRE2_pile2_oi2
#   
#   pair_id <- pair_id+1
#   
# }
# 
# 
# df_mBC_w_CRE2_pile2_conn_clu <- rbindlist(df_list)
# 
# df_mBC_w_CRE2_pile_max_conn <- df_mBC_w_CRE2_pile2_conn_clu %>% group_by(conn_clu_id) %>% slice_max(mBC_counts) 
# 
# # ggplot(df_mBC_w_CRE2_pile) + stat_bin(aes(x=mBC_counts),geom="step")+scale_x_log10()+scale_y_log10()
# # ggplot(df_mBC_w_CRE2_pile_max_conn) + stat_bin(aes(x=mBC_counts),geom="step")+scale_x_log10()+scale_y_log10()
# 
# mBC_max_read_thresh <- 2
# ggplot() + stat_bin(aes(x=df_mBC_w_CRE2_pile_max_conn$mBC_counts),geom="step",color="red")+
#   stat_bin(aes(x=df_mBC_w_CRE2_pile$mBC_counts),geom="step",color="black")+scale_x_log10()+scale_y_log10()+
#   geom_vline(aes(xintercept=mBC_max_read_thresh),linetype="dotted")
# 
# 
# # final threshold!
# df_mBC_w_CRE2_pile_max_conn_hi <- df_mBC_w_CRE2_pile_max_conn %>% filter(mBC_counts>=mBC_max_read_thresh)
# dim(df_mBC_w_CRE2_pile_max_conn_hi)
# 

# 
# 
# # manual remove unexpected CREs
# # df_mBC_w_CRE2_pile_max_conn_hi <- df_mBC_w_CRE2_pile_max_conn_hi %>% 
# #   filter(!(CRE_R_id %in% c("Gata4_chr14_5729","Epas1_chr17_10063","Col5a1_chr2_2586"))) %>%
# #   filter(!(CRE_L_id %in% c("Cdk5r1_chr11_12590","Sox2_chr3_2007","Sox2_chr3_2009")))
# 
# df_mBC_w_CRE2_pile_max_conn_hi <- df_mBC_w_CRE2_pile_max_conn_hi %>% 
#   filter(!(CRE_L_id %in% c("Gata4_chr14_5729","Epas1_chr17_10063","Col5a1_chr2_2586"))) %>%
#   filter(!(CRE_R_id %in% c("Cdk5r1_chr11_12590","Sox2_chr3_2007","Sox2_chr3_2009")))
# 
# final_out <- sprintf("final_suba_%s_%s.txt",output_file_id,date_str)
# write.table(df_mBC_w_CRE2_pile_max_conn_hi,final_out,
#             sep="\t",quote=FALSE,row.names=FALSE)
# 
# table(df_mBC_w_CRE2_pile_max_conn_hi$CRE_L_id,
#       df_mBC_w_CRE2_pile_max_conn_hi$CRE_R_id)
# 
# table(paste0(df_mBC_w_CRE2_pile_max_conn_hi$CRE_L_id,df_mBC_w_CRE2_pile_max_conn_hi$CRE_L_ori),
#       paste0(df_mBC_w_CRE2_pile_max_conn_hi$CRE_R_id,df_mBC_w_CRE2_pile_max_conn_hi$CRE_R_ori))
# 
# 
# # QC metrics/numbers:
# df_aln2 %>% pull(target_name) %>% unique %>% length
# sum(nanop_read_l>5000 & nanop_read_l<10000)
# df_aln2_hi2 %>% pull(target_name) %>% unique %>% length
# dim(df_usable_reads)[1]
# dim(df_mBC_w_CRE2_pile2)[1]
# dim(df_mBC_w_CRE2_pile_max_conn_hi)[1]
# 
# 
# ## singling specific reads as example:
# 
# read_ids_example_CRE_pairs <- 
#   df_mBC_w_CRE2 %>% filter(CRE_L_id=="Sox2_chr3_2007" & 
#                              CRE_L_ori=="-" & 
#                              CRE_R_id=="Epas1_chr17_10063" & 
#                              CRE_R_ori=="-") %>% pull(target_name)
# 
# example_reads_pairs <- df_fastq[df_fastq@id %in% read_ids_example_CRE_pairs[1:100]]@sread
# 
# write.fasta(as.list(example_reads_pairs),names=read_ids_example_CRE_pairs[1:100],
#             sprintf("%s_example_Sox2_chr3_2007-_Epas1_chr17_10063-_reads_BC1_%s.fasta",output_file_id,date_str))
