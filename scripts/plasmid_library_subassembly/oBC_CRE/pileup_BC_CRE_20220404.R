
library(dplyr)
library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
out_file_pileup <- args[2]
out_file_full_pileup <- args[3]


# read in BC-CRE merged table
df_BC_CRE <- read.table(input_file, sep="\t",header=TRUE)


# pile up, keeping unique tagmentation position information, thresholding on mapping quality
mapq_thresh <- 40

df_BC_CRE_pileup <- df_BC_CRE %>% filter(mapq>mapq_thresh) %>% group_by(BC,CRE_id,read_forward_in_ref,read_start) %>% summarize(n_counts=length(read_start)) 
df_BC_CRE_pileup2 <- df_BC_CRE_pileup %>% arrange(CRE_id,BC)
df_BC_CRE_pileup3 <- df_BC_CRE_pileup2 %>% transform(BC_CRE_strand=paste0(BC,"_",CRE_id,"_",read_forward_in_ref))

write.table(df_BC_CRE_pileup3,out_file_pileup,
            sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

# pile up also over position information. calculate summary stats of position dispersion as a downstream flag 
df_BC_CRE_full_pileup <- df_BC_CRE %>% filter(mapq>mapq_thresh) %>% group_by(BC,CRE_id,read_forward_in_ref) %>% 
  summarize(n_diff_pos=length(unique(read_start)), 
            n_total_counts=length(read_start),
            q10_pos=quantile(read_start,0.1),
            q50_pos=quantile(read_start,0.5),
            q90_pos=quantile(read_start,0.9)) %>%  arrange(CRE_id,BC) %>% 
  transform(BC_CRE_strand=paste0(BC,"_",CRE_id,"_",read_forward_in_ref))

df_counts_per_BC <- df_BC_CRE_full_pileup %>% group_by(BC) %>% summarize(tot_reads_per_BC=sum(n_total_counts))
df_BC_CRE_full_pileup2 <- df_BC_CRE_full_pileup %>% left_join(df_counts_per_BC) %>% transform(prop_reads_per_CRE=n_total_counts/tot_reads_per_BC)


write.table(df_BC_CRE_full_pileup2,out_file_full_pileup,
            sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
