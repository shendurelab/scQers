
library(dplyr)
library(tidyverse)
library(scales)
library(Biostrings)
library(stringr)
library(patchwork)


# input arguments
args <- commandArgs(trailingOnly=TRUE)
file_oi <- args[1] # input file (full pileup from pileup_BC_CRE_20220404.R output)
sample_oi <- args[2]  # file name
read_thresh <- args[3] # threshold read count
file_oBC_mBC_pair <- args[4]  # previously associated oBC-mBC pairs to complete triplet
QC_plot_bool <- args[5]   # plot QC or no? 


setwd("/Users/jbl/Documents/UW/manuscripts/single_cell_expression_reporter/addgene_lib_submission/BC_CRE_processing")
file_oi <- "example_CRE_oBC_full_pileup.txt.gz"
sample_oi <- "example_CRE_oBC"
read_thresh <- 3
file_oBC_mBC_pair <- "p025_recloned_complex_mBC_oBC_subassembly.txt.gz"
QC_plot_bool <- TRUE

# a priori two barcode pairs (previously obtained)
df_mBC_oBC_pairs <- read.table(file_oBC_mBC_pair,sep="\t",header=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # 
# thresholds for selecting valid associations
# # # # # # # # # # # # # # # # # # # # # # # # 

prop_thresh <- 0.975  # proportion of reads from a BC associated to a given CRE

# positional threshold and variability: related to cut size of the semi-specific PCR product
median_pos_thresh_lo <- 30
median_pos_thresh_hi <- 400   
pos_90_to_10_thresh_lo <- 30
pos_90_to_10_thresh_hi <- 500



# read full pile up info
df_full_pile_up3 <- read.table(file_oi, sep="\t", header=TRUE) %>% 
  filter(!str_detect(CRE_id,"constant_backbone")) 

# repeat proportion and count calculation on the variant without the backbone pieces:
df_full_pile_up_high <- df_full_pile_up3 %>% 
  filter(n_total_counts>read_thresh & read_forward_in_ref=="True")


# filtering by expected tagmentation positions within CRE & orientation
df_full_pile_up_pos_cut <- df_full_pile_up3 %>% filter(q50_pos<median_pos_thresh_hi &
                                                         q50_pos>median_pos_thresh_lo & 
                                                         read_forward_in_ref=="True")

df_counts_per_BC <- df_full_pile_up_pos_cut %>% group_by(BC) %>% 
  summarize(tot_reads_per_BC=sum(n_total_counts))

df_full_pile_up_pos_cut2 <- df_full_pile_up_pos_cut %>% 
  left_join(df_counts_per_BC) %>% 
  transform(prop_reads_per_CRE=n_total_counts/tot_reads_per_BC)


df_full_pile_up_high2 <- df_full_pile_up_pos_cut2 %>% filter( 
  (q90_pos-q10_pos)<median_pos_thresh_hi & 
    (q90_pos-q10_pos)>median_pos_thresh_lo & 
    n_total_counts>read_thresh &
    prop_reads_per_CRE>prop_thresh)





# join w/ previously determimed oBC-mBC pairing
df_full_pile_up_high_w_mBC <- df_full_pile_up_high2 %>% 
  rename(c("BC"="RC_oBC", "BC_CRE_strand"="RC_oBC_CRE_strand","tot_reads_per_BC"="tot_reads_per_oBC"))

RC_oBC <- DNAStringSet(df_full_pile_up_high_w_mBC$RC_oBC)
oBC <- reverseComplement(RC_oBC)
df_full_pile_up_high_w_mBC$oBC <- as.character(oBC)
df_full_pile_up_high_w_mBC <- df_full_pile_up_high_w_mBC %>% left_join(df_mBC_oBC_pairs)


# summarizing per CRE
df_CRE_summary <- df_full_pile_up_high2 %>% group_by(CRE_id) %>% summarize(n_reads_per_CRE=sum(n_total_counts), 
                                                                           n_BC_per_CRE=length(unique(BC)))

write.table(df_full_pile_up_high_w_mBC,
            sprintf("%s_final_triplets.txt",sample_oi),
            sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(df_CRE_summary,
            sprintf("%s_summary_per_CRE.txt",sample_oi),
            sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)



if (QC_plot_bool){
  # # # # # # # # # # # # 
  # QC plots
  # # # # # # # # # # # # 
  
  labels_reads <- 10^(0:4)
  lims_reads <- c(1,10^4)
  labels_pair_counts <- 10^(0:5)
  lims_pair_counts <- c(1,10^5)
  plt_read_dist <- ggplot(df_full_pile_up3)+stat_bin(aes(x=n_total_counts),geom="step",bins=50)+
    geom_vline(xintercept=read_thresh,color="red",linetype="dashed")+
    scale_x_log10(breaks = labels_reads,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_reads)+
    scale_y_log10(breaks = labels_pair_counts,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_pair_counts)+annotation_logticks()+
    labs(x="Reads per CRE-oBC pair",y="Count (PDF)")
  # plt_read_dist
  
  labels_tag_pos <- 10^(0:3)
  lims_tag_pos <- c(1,500)
  plt_tag_pos_dist <- ggplot(df_full_pile_up3)+stat_bin(aes(x=n_diff_pos),geom="step",bins=50)+
    scale_x_log10(breaks = labels_tag_pos,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_tag_pos)+
    scale_y_log10(breaks = labels_pair_counts,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_pair_counts)+annotation_logticks()+
    labs(x="Distinct number of tagmented positions\n per CRE-oBC pair (~UMI)",y="Count (PDF)")
  # plt_tag_pos_dist
  
  plt_reads_vs_tag_pos<- ggplot(df_full_pile_up3)+geom_hex(aes(x=n_total_counts,y=n_diff_pos,fill=log(..count..)))+
    scale_y_log10(breaks = labels_tag_pos,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_tag_pos)+
    scale_x_log10(breaks = labels_reads,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_reads)+annotation_logticks()+
    geom_abline(color="grey",linetype="dashed")+
    labs(y="Distinct number of tagmented positions\n per CRE-oBC pair (~UMI)",x="Reads per CRE-oBC pair")
  # plt_reads_vs_tag_pos
  
  
  plt_tag_pos_high_reads <- ggplot(df_full_pile_up_high)+geom_hex(aes(x=q50_pos,y=q90_pos-q10_pos))+
    geom_vline(xintercept=c(median_pos_thresh_lo,median_pos_thresh_hi),color="red",linetype="dashed")+
    geom_hline(yintercept=c(pos_90_to_10_thresh_lo,pos_90_to_10_thresh_hi),color="red",linetype="dashed")+
    labs(x="Median tagmented positions within CRE",y="Spread tagmented position (10-to-90th percentile) ")
  # plt_tag_pos_high_reads
  
  plt_reads_vs_prop <- ggplot(df_full_pile_up_pos_cut2)+geom_hex(aes(x=tot_reads_per_BC,y=prop_reads_per_CRE,fill=log(..count..)),bins=100)+ 
    scale_x_log10(breaks = labels_reads,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims_reads)+annotation_logticks(side="b")+
    geom_hline(yintercept=prop_thresh,color="red",linetype="dashed")+geom_vline(xintercept=read_thresh,color="red",linetype="dashed")+
    labs(y="Proportion of reads from oBC to CRE",x="Reads per CRE-oBC pair")
  # plt_reads_vs_prop
  
  
  lim_reads_per_CRE <- c(1E3,5E5)
  labels_BC_per_CRE <- c(1000,1E4,1E5)
  
  plt_tot_reads_per_CRE <- ggplot(df_CRE_summary) + stat_bin(aes(x=n_reads_per_CRE),geom="step",bins=25)+
    scale_x_log10(breaks = labels_BC_per_CRE,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lim_reads_per_CRE )+annotation_logticks(side="b")+
    labs(y="Coutn (PDF)",x="Reads per CRE \n(across valid oBC-CRE pairs)")
  plt_tot_reads_per_CRE
  
  lim_BC_per_CRE <- c(1,1E4)
  labels_BC_per_CRE <- c(1,10,100,1000)
  
  plt_num_final_BC_per_CRE <- ggplot(df_CRE_summary) + stat_bin(aes(x=n_BC_per_CRE),geom="step")+
    scale_x_log10(breaks = labels_BC_per_CRE,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lim_BC_per_CRE)+annotation_logticks(side="b")+
    labs(y="Count (PDF)",x="oBC per CRE \n(across valid oBC-CRE pairs)",
         title=sprintf('n=%d valid pairs\n',dim(df_full_pile_up_high_w_mBC)[1]))
  plt_num_final_BC_per_CRE
  
  lims_pair_counts2 <- c(1,1E4)
  labels_pair_counts2 <- 10^(0:4)
  
  plt_dist_prop_high_read_pairs <- ggplot(df_full_pile_up_pos_cut2 %>% filter(n_total_counts>read_thresh))+stat_bin(aes(y=prop_reads_per_CRE),bins=100)+ 
    geom_hline(yintercept=prop_thresh,color="red",linetype="dashed")+ylim(c(0,1.02))+
    scale_x_log10(breaks = labels_pair_counts2,
                  labels = trans_format("log10", math_format(10^.x)))+annotation_logticks(side="b")+
    labs(y="Proportion of reads from oBC to CRE",x="Number oBC-CRE pairs\n (above read threshold)")
  plt_dist_prop_high_read_pairs
  
  final_plt <- plt_read_dist+
    plt_tag_pos_dist+
    plt_reads_vs_tag_pos+
    plot_spacer()+
    plot_spacer()+
    plot_spacer()+
    plt_reads_vs_prop+
    plt_dist_prop_high_read_pairs+
    plt_tag_pos_high_reads+
    plot_spacer()+
    plot_spacer()+
    plot_spacer()+
    plt_tot_reads_per_CRE+
    plt_num_final_BC_per_CRE+
    plot_layout(ncol=3,
                heights = c(4,1,4,1,4)) & theme(title = element_text(size = 6),
                                                axis.text.x = element_text(size=6),
                                                axis.text.y = element_text(size=6),
                                                legend.text = element_text(size=6),
                                                legend.title = element_text(size=6))
  # dev.new()
  # final_plt
  
  
  
  fig_name <- sprintf("%s_summary_plot.pdf",sample_oi)
  pdf(fig_name,height=10,width=10)
  print(final_plt)
  dev.off()
}







