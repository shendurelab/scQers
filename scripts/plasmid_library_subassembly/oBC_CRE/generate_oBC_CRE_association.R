
library(dplyr)
library(tidyverse)
library(scales)
library(Biostrings)
library(stringr)
library(Seurat)
library(Matrix)
library(reshape2)
library(patchwork)
library(igraph)
library(pracma)



dir_outs <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq023_devCRE_subassemb_sci_mBC_oBC_MPRA_20220331/oBC_devCRE_subassembly_final_outs_w_multimappers_20220818/"
setwd(dir_outs)


# load full list of CREs: 
df_all_loci <- read.table("/Users/jbl/Dropbox (MIT)/sc_EB_RNA_ATAC_plots/devCRE_PCR_cloning_analysis_20211031/devCRE_region_seq_no_seq_primers_v2_20220209.txt",
                          sep="\t",header=TRUE)
df_all_loci <- df_all_loci %>% select(CRE_id=peak_id, everything())


# p25 1x
file_oBC_mBC_pair <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq020_p25_asso_p22_mBC_sc_rep_CS2_v2_bulkoBC_sci_enrich_20220321/p25_mBC_oBC/final_table_p25_v2_oBC_mBC_subassembly_p25_mBC_oBC_asso_1x_S14_20220320.txt"


sample_oi <- "p55_1x_10m2_oBC_devCRE_suba_S2_w_multimapper"
file_oi <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq023_devCRE_subassemb_sci_mBC_oBC_MPRA_20220331/pile_up_oBC_devCRE_overlap_separated/p55_1x_10m2_oBC_devCRE_suba_S2_full_pileup_k2_20220818.txt.gz"
read_thresh <- 30
prop_thresh <- 0.95

median_pos_thresh_hi <- 300
pos_90_to_10_thresh_hi <- 300
median_pos_thresh_lo <- 30
pos_90_to_10_thresh_lo <- 30



df_mBC_oBC_pairs <- read.table(file_oBC_mBC_pair,sep="\t",header=TRUE)


# read full pile up info
df_full_pile_up3 <- read.table(file_oi, sep="\t", header=TRUE)



# # QC plots
# labels_reads <- 10^(0:4)
# lims_reads <- c(1,10^4)
# labels_pair_counts <- 10^(0:5)
# lims_pair_counts <- c(1,10^5)
# plt_read_dist <- ggplot(df_full_pile_up3)+stat_bin(aes(x=n_total_counts),geom="step",bins=50)+
#   geom_vline(xintercept=read_thresh,color="red",linetype="dashed")+
#   scale_x_log10(breaks = labels_reads,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_reads)+
#   scale_y_log10(breaks = labels_pair_counts,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_pair_counts)+annotation_logticks()+
#   labs(x="Reads per CRE-oBC pair",y="Count (PDF)")
# plt_read_dist
# 
# labels_tag_pos <- 10^(0:3)
# lims_tag_pos <- c(1,500)
# plt_tag_pos_dist <- ggplot(df_full_pile_up3)+stat_bin(aes(x=n_diff_pos),geom="step",bins=50)+
#   scale_x_log10(breaks = labels_tag_pos,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_tag_pos)+
#   scale_y_log10(breaks = labels_pair_counts,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_pair_counts)+annotation_logticks()+
#   labs(x="Distinct number of tagmented positions\n per CRE-oBC pair (~UMI)",y="Count (PDF)")
# plt_tag_pos_dist
# 
# 
# plt_reads_vs_tag_pos<- ggplot(df_full_pile_up3)+geom_hex(aes(x=n_total_counts,y=n_diff_pos,fill=log(..count..)))+
#   scale_y_log10(breaks = labels_tag_pos,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_tag_pos)+
#   scale_x_log10(breaks = labels_reads,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_reads)+annotation_logticks()+
#   geom_abline(color="grey",linetype="dashed")+
#   labs(y="Distinct number of tagmented positions\n per CRE-oBC pair (~UMI)",x="Reads per CRE-oBC pair")
# plt_reads_vs_tag_pos






#  proportion and count calculation on the variant without the backbone pieces:

df_full_pile_up_high <- df_full_pile_up3 %>% filter(n_total_counts>read_thresh & read_forward_in_ref=="True")
plt_tag_pos_high_reads <- ggplot(df_full_pile_up_high)+geom_hex(aes(x=q50_pos,y=q90_pos-q10_pos))+
  geom_vline(xintercept=c(median_pos_thresh_lo,median_pos_thresh_hi),color="red",linetype="dashed")+
  geom_hline(yintercept=c(pos_90_to_10_thresh_lo,pos_90_to_10_thresh_hi),color="red",linetype="dashed")+
  labs(x="Median tagmented positions within CRE",y="Spread tagmented position (10-to-90th percentile) ")
plt_tag_pos_high_reads


# adding positional cutoffs
df_full_pile_up_pos_cut <- df_full_pile_up3 %>% filter(q50_pos<median_pos_thresh_hi &
                                                         q50_pos>median_pos_thresh_lo & 
                                                         read_forward_in_ref=="True")

                                                    
                                                    
df_full_pile_up_pos_cut2 <- df_full_pile_up_pos_cut %>% filter(!str_detect(CRE_id,"constant_backbone")) %>% select(-c("tot_reads_per_BC","prop_reads_per_CRE"))
df_counts_per_BC <- df_full_pile_up_pos_cut2 %>% group_by(BC) %>% summarize(tot_reads_per_BC=sum(n_total_counts))
df_full_pile_up_pos_cut3 <- df_full_pile_up_pos_cut2 %>% left_join(df_counts_per_BC) %>% transform(prop_reads_per_CRE=n_total_counts/tot_reads_per_BC)



# # QC plots
# plt_reads_vs_prop <- ggplot(df_full_pile_up_pos_cut3)+geom_hex(aes(x=tot_reads_per_BC,y=prop_reads_per_CRE,fill=log(..count..)),bins=100)+ 
#   scale_x_log10(breaks = labels_reads,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lims_reads)+annotation_logticks(side="b")+
#   geom_hline(yintercept=prop_thresh,color="red",linetype="dashed")+geom_vline(xintercept=read_thresh,color="red",linetype="dashed")+
#   labs(y="Proportion of reads from oBC to CRE",x="Reads per CRE-oBC pair")
# plt_reads_vs_prop


df_full_pile_up_high2 <- df_full_pile_up_pos_cut3 %>% filter( 
  (q90_pos-q10_pos)<median_pos_thresh_hi & 
    (q90_pos-q10_pos)>median_pos_thresh_lo & 
    n_total_counts>read_thresh &
    prop_reads_per_CRE>prop_thresh & 
    !(CRE_id %in% short_hom_CREs))

# special case with short homology to each other breaking the call from the small proportion mapping to each other.
# seems to be from a repetitive element as that sequence is all over the place in the genome. 
short_hom_CREs <- c("Gata4_chr14_5749", "Txndc12_chr4_7975")

# first one
df_full_pile_up_pos_cut2_oi1 <- df_full_pile_up_pos_cut %>% 
  filter(!(str_detect(CRE_id,"constant_backbone") | CRE_id==short_hom_CREs[1])) %>%
  select(-c("tot_reads_per_BC","prop_reads_per_CRE"))
df_counts_per_BC_oi1 <- df_full_pile_up_pos_cut2_oi1 %>% group_by(BC) %>% summarize(tot_reads_per_BC=sum(n_total_counts))
df_full_pile_up_pos_cut3_oi1 <- df_full_pile_up_pos_cut2_oi1 %>% left_join(df_counts_per_BC_oi1) %>% transform(prop_reads_per_CRE=n_total_counts/tot_reads_per_BC)

df_full_pile_up_high_oi1 <- df_full_pile_up_pos_cut3_oi1 %>% filter( 
  (q90_pos-q10_pos)<median_pos_thresh_hi & 
    (q90_pos-q10_pos)>median_pos_thresh_lo & 
    n_total_counts>read_thresh &
    prop_reads_per_CRE>prop_thresh & 
    CRE_id==short_hom_CREs[2])

df_full_pile_up_pos_cut2_oi2 <- df_full_pile_up_pos_cut %>% 
  filter(!(str_detect(CRE_id,"constant_backbone") | CRE_id==short_hom_CREs[2])) %>%
  select(-c("tot_reads_per_BC","prop_reads_per_CRE"))
df_counts_per_BC_oi2 <- df_full_pile_up_pos_cut2_oi2 %>% group_by(BC) %>% summarize(tot_reads_per_BC=sum(n_total_counts))
df_full_pile_up_pos_cut3_oi2 <- df_full_pile_up_pos_cut2_oi2 %>% left_join(df_counts_per_BC_oi2) %>% transform(prop_reads_per_CRE=n_total_counts/tot_reads_per_BC)

df_full_pile_up_high_oi2 <- df_full_pile_up_pos_cut3_oi2 %>% filter( 
  (q90_pos-q10_pos)<median_pos_thresh_hi & 
    (q90_pos-q10_pos)>median_pos_thresh_lo & 
    n_total_counts>read_thresh &
    prop_reads_per_CRE>prop_thresh & 
    CRE_id==short_hom_CREs[1])


df_full_pile_up_high3 <- rbind(df_full_pile_up_high2,df_full_pile_up_high_oi1,df_full_pile_up_high_oi2)



# # QC plots
# lims_pair_counts2 <- c(1,1E4)
# labels_pair_counts2 <- 10^(0:4)
# plt_dist_prop_high_read_pairs <- ggplot(df_full_pile_up_pos_cut3 %>% filter(n_total_counts>read_thresh))+stat_bin(aes(y=prop_reads_per_CRE),bins=100)+ 
#   geom_hline(yintercept=prop_thresh,color="red",linetype="dashed")+ylim(c(0,1.02))+
#   scale_x_log10(breaks = labels_pair_counts2,
#                 labels = trans_format("log10", math_format(10^.x)))+annotation_logticks(side="b")+
#   labs(y="Proportion of reads from oBC to CRE",x="Number oBC-CRE pairs\n (above read threshold)")
# plt_dist_prop_high_read_pairs



# join w/ oBC-mBC pairing
df_full_pile_up_high_w_mBC <- df_full_pile_up_high3 %>% rename(c("BC"="RC_oBC", "BC_CRE_strand"="RC_oBC_CRE_strand","tot_reads_per_BC"="tot_reads_per_oBC"))

RC_oBC <- DNAStringSet(df_full_pile_up_high_w_mBC$RC_oBC)
oBC <- reverseComplement(RC_oBC)
df_full_pile_up_high_w_mBC$oBC <- as.character(oBC)
df_full_pile_up_high_w_mBC <- df_full_pile_up_high_w_mBC %>% left_join(df_mBC_oBC_pairs)


# summarizing per CRE
df_CRE_summary <- df_full_pile_up_high3 %>% group_by(CRE_id) %>% summarize(n_reads_per_CRE=sum(n_total_counts), 
                                                                                 n_BC_per_CRE=length(unique(BC)))

df_all_loci_w_mapping_info <- df_all_loci %>% left_join(df_CRE_summary) %>% 
  mutate(n_reads_per_CRE=replace_na(n_reads_per_CRE,0)) %>% 
  mutate(n_BC_per_CRE=replace_na(n_BC_per_CRE,0))


# # QC plots
# lim_reads_per_CRE <- c(1E3,5E5)
# labels_BC_per_CRE <- c(1000,1E4,1E5)
# 
# plt_tot_reads_per_CRE <- ggplot(df_all_loci_w_mapping_info) + stat_bin(aes(x=n_reads_per_CRE+50),geom="step",bins=25)+
#   scale_x_log10(breaks = labels_BC_per_CRE,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lim_reads_per_CRE )+annotation_logticks(side="b")+
#   labs(y="Coutn (PDF)",x="Reads per CRE \n(across valid oBC-CRE pairs)")
# plt_tot_reads_per_CRE
# 
# lim_BC_per_CRE <- c(1,1E4)
# labels_BC_per_CRE <- c(1,10,100,1000)
# 
# plt_num_final_BC_per_CRE <- ggplot(df_all_loci_w_mapping_info) + stat_bin(aes(x=n_BC_per_CRE),geom="step")+
#   scale_x_log10(breaks = labels_BC_per_CRE,
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = lim_BC_per_CRE)+annotation_logticks(side="b")+
#   labs(y="Count (PDF)",x="oBC per CRE \n(across valid oBC-CRE pairs)",
#        title=sprintf('n=%d valid pairs\nn=%d complete dropouts',dim(df_full_pile_up_high_w_mBC)[1],
#                      sum(df_all_loci_w_mapping_info$n_reads_per_CRE==0)))
# plt_num_final_BC_per_CRE
# 
# 
# final_plt <- plt_read_dist+
#   plt_tag_pos_dist+
#   plt_reads_vs_tag_pos+
#   plot_spacer()+
#   plot_spacer()+
#   plot_spacer()+
#   plt_reads_vs_prop+
#   plt_dist_prop_high_read_pairs+
#   plt_tag_pos_high_reads+
#   plot_spacer()+
#   plot_spacer()+
#   plot_spacer()+
#   plt_tot_reads_per_CRE+
#   plt_num_final_BC_per_CRE+
#   plot_layout(ncol=3,
#               heights = c(4,1,4,1,4)) & theme(title = element_text(size = 6),
#                               axis.text.x = element_text(size=6),
#                               axis.text.y = element_text(size=6),
#                               legend.text = element_text(size=6),
#                               legend.title = element_text(size=6))
# 
# 
# fig_name <- sprintf("%s_summary_plot_20220818.pdf",sample_oi)
# pdf(fig_name,height=10,width=10)
# print(final_plt)
# dev.off()

valid_reads <- sum(df_full_pile_up_high_w_mBC$n_total_counts)
all_reads <- sum(df_full_pile_up3$n_total_counts)
df_meta_summary_oi <- data.frame(sample=sample_oi,
                                 read_thresh,prop_thresh,
                                 n_valid_oBC_CRE=dim(df_full_pile_up_high_w_mBC)[1],
                                 median_valid_oBC_per_CRE=median(df_all_loci_w_mapping_info$n_BC_per_CRE),
                                 number_CRE_dropouts=sum(df_all_loci_w_mapping_info$n_BC_per_CRE==0),
                                 fraction_aligned_reads_to_valid_construct=valid_reads/all_reads)
df_meta_summary <- rbind(df_meta_summary,df_meta_summary_oi)




#  output files 

write.table(df_full_pile_up_high_w_mBC,
            sprintf("%s_final_triplets_20220818.txt",sample_oi),
            sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(df_all_loci_w_mapping_info,
            sprintf("%s_summary_per_CRE_20220818.txt",sample_oi),
            sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)


