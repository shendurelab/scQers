
library(tidyverse)
library(stringr)
library(ArchR)
library(ggnewscale)
library(cowplot)
library(BSgenome)

# loading annotation
regions <- readRDS("dependencies/mm10_gene_anno_ArchR.RDS")

# output from TFBS_locus_scan
path <- "locus_wide_TFBS_all_loci_20220918"
files <- dir(path)

df_devCREs <- read.table("dependencies/target_devCRE_peak_set_20210803.txt", 
                         sep="\t",header=TRUE)


genes <- str_replace(str_replace(files,"TFBS_hits_",""),"_100kb_pm_TSS_20220915.txt","")
names(files) <- genes



size_region <- 2E5+1
df_pos <- data.frame(start_pos=seq(1,size_region))
averaging_size <- 500
mask_dist <- 500



affinity_cutoff <- 0.4

df_kmer_w_hits_masked_windows_wd_all <- data.frame()
hits_TFBS_devCRE_all <- data.frame()
for (gene_oi in genes){
  
  print(gene_oi)
  
  
  df_TFBS_hits <- read.table(sprintf("%s/%s",path,files[gene_oi]),header=TRUE)
  df_TFBS_hits <- df_TFBS_hits %>% filter(value>=affinity_cutoff)
  
  
  # annotating all positions with TFBS hit positions, with a binarized TFBS hit score
  df_kmer_w_hits <- data.frame()
  TF_names <- c("Gata4","Foxa2","Sox17")
  for (TF_oi in TF_names){
    df_kmer_w_hits_oi <- df_pos
    df_kmer_w_hits_oi$TF_name <- TF_oi
    df_kmer_w_hits_oi$hits <- FALSE
    df_kmer_w_hits_oi$affinity <- NA
    df_kmer_w_hits_oi$hits[df_TFBS_hits %>% filter(TF_name==TF_oi) %>% pull(start_pos)] <- TRUE
    df_kmer_w_hits_oi$affinity[df_TFBS_hits %>% filter(TF_name==TF_oi) %>% pull(start_pos)] <- df_TFBS_hits %>% filter(TF_name==TF_oi) %>% pull(value)
    
    df_kmer_w_hits <- rbind(df_kmer_w_hits,df_kmer_w_hits_oi)
  }
  
  df_kmer_w_hits2 <- df_kmer_w_hits
  df_kmer_w_hits2$binarized_score <- 0
  df_kmer_w_hits2$binarized_score[df_kmer_w_hits2$hits] <- 1
  
  df_kmer_w_hits2$number_TFBS_window <- NA
  
  for (TF_oi in TF_names){

    thresholded_score <- df_kmer_w_hits2 %>% filter(TF_name==TF_oi) %>% pull(binarized_score)
    sum_hits_region <- stats:::filter(thresholded_score,rep(1,averaging_size),method = "convolution", sides=1)
    df_kmer_w_hits2$number_TFBS_window[df_kmer_w_hits2$TF_name==TF_oi] <- sum_hits_region
  }
  
  
  # extract counts in the CREs from locus of interest
  df_devCREs_oi <- df_devCREs %>% 
    filter(DE_gene==gene_oi) %>% 
    transform(CRE_id = sprintf("%s_%s",gene_oi,unique_peak_id), mid_point=0.5*(start+end)) %>% 
    select(CRE_id, genome_positions=start,mid_point)
  
  offset_pos <- ceiling(averaging_size/2)
  df_devCREs_oi$genome_positions <- df_devCREs_oi$mid_point+offset_pos
  
  
  # masking around CREs to avoid biasing results, +/- averaging window around midpoint
  # position info for gene of interest
  region_oi <- regions[which(tolower(mcols(regions)$symbol) %in% tolower(gene_oi))]
  region_oi <- region_oi[order(match(tolower(mcols(region_oi)$symbol), tolower(gene_oi)))]
  region_oi_TSS <- resize(region_oi, 1, "start")
  strand(region_oi_TSS) <- "*"
  peak_oi_range_ext_long2 <- extendGR(region_oi_TSS, upstream = ext_size_long2, downstream = ext_size_long2)
  
  chr_oi <- seqnames(peak_oi_range_ext_long2) %>% as.character()
  pos_oi <- start(peak_oi_range_ext_long2):end(peak_oi_range_ext_long2)
  pos_df <- data.frame(chr=chr_oi,genome_positions=pos_oi,relative_positions=seq(length(pos_oi)))
  
  inds_masked <- c()
  for (mid_CRE_pos in df_devCREs_oi$mid_point){
    inds_masked <- c(inds_masked, (mid_CRE_pos-mask_dist):(mid_CRE_pos+mask_dist))
  }
  
  
  df_kmer_w_hits3 <- df_kmer_w_hits2 %>% left_join(pos_df %>% select(chr,genome_positions,start_pos=relative_positions))
  df_kmer_w_hits_masked <- df_kmer_w_hits3
  df_kmer_w_hits_masked2 <- df_kmer_w_hits_masked %>% left_join(pos_df %>% select(chr,genome_positions,start_pos=relative_positions))
  df_kmer_w_hits_masked2$number_TFBS_window_CRE_masked <- df_kmer_w_hits_masked2$number_TFBS_window
  df_kmer_w_hits_masked2$number_TFBS_window_CRE_masked[df_kmer_w_hits_masked2$genome_positions %in% inds_masked] <- NA
  
  
  # only get data in windows (steps of half the averaging window size)
  df_kmer_w_hits_masked_windows <- df_kmer_w_hits_masked2 %>% filter( (start_pos %% (averaging_size/2))==0) %>% select(-c("binarized_score","affinity","number_TFBS_window","hits"))
  df_kmer_w_hits_masked_windows_wd <- df_kmer_w_hits_masked_windows %>% pivot_wider(names_from=TF_name,values_from=number_TFBS_window_CRE_masked)
  
  
  hits_TFBS_devCRE <- df_devCREs_oi %>% 
    left_join(df_kmer_w_hits3 %>% select(-c("binarized_score","affinity","hits","start_pos"))) %>% 
    pivot_wider(names_from=TF_name,values_from=number_TFBS_window)
  
  
  df_kmer_w_hits_masked_windows_wd_all <- rbind(df_kmer_w_hits_masked_windows_wd_all,
                                                df_kmer_w_hits_masked_windows_wd %>% transform(locus=gene_oi))
  
  hits_TFBS_devCRE_all <- rbind(hits_TFBS_devCRE_all,
                                hits_TFBS_devCRE %>% transform(locus=gene_oi))
  
  
}



# layer information about distal/proximal/specific:
df_sc_meta_data <- read.table("dependencies/scQer_activity_specificity_aggregated_lvl2_final_list_w_metadata_20220905.txt",header=TRUE)
hits_TFBS_devCRE_all2 <- hits_TFBS_devCRE_all %>% left_join(df_sc_meta_data %>% select(CRE_id,is_promoter,all_rep_spec_hit,all_rep_act_hit))


# color scheme
hits_TFBS_devCRE_all2$category2 <- "inactive"
hits_TFBS_devCRE_all2$category2[hits_TFBS_devCRE_all2$all_rep_act_hit & !hits_TFBS_devCRE_all2$all_rep_spec_hit & hits_TFBS_devCRE_all2$is_promoter] <- "non-specific, active, promoter"
hits_TFBS_devCRE_all2$category2[hits_TFBS_devCRE_all2$all_rep_act_hit & !hits_TFBS_devCRE_all2$all_rep_spec_hit & !hits_TFBS_devCRE_all2$is_promoter] <- "non-specific, active"
hits_TFBS_devCRE_all2$category2[hits_TFBS_devCRE_all2$all_rep_spec_hit] <- "specific, active"
table(hits_TFBS_devCRE_all2$category2)
col_dict2 <- c("inactive"="grey90","non-specific, active, promoter" = "goldenrod1","non-specific, active"="black","specific, active"="red")



df_quantiles <- df_kmer_w_hits_masked_windows_wd_all %>% group_by(locus) %>% 
  summarize(q90_Gata4=quantile(Gata4,0.9,na.rm=TRUE),
            q95_Gata4=quantile(Gata4,0.95,na.rm=TRUE),
            q99_Gata4=quantile(Gata4,0.99,na.rm=TRUE),
            max_Gata4=max(Gata4,na.rm=TRUE),
            q90_Sox17=quantile(Sox17,0.9,na.rm=TRUE),
            q95_Sox17=quantile(Sox17,0.95,na.rm=TRUE),
            q99_Sox17=quantile(Sox17,0.99,na.rm=TRUE),
            max_Sox17=max(Sox17,na.rm=TRUE)) %>% filter(locus %in% endo_loci)



endo_loci <- c("Bend5","Epas1","Foxa2","Gata4","Lama1","Lamb1","Lamc1","Sparc","Btg1","Cited2","Klf4","Sox17","Txndc12")


df_locus_order <- data.frame(locus=c(endo_loci),
                             locus_order=seq(1:13))

df_kmer_w_hits_masked_windows_wd_all2 <- df_kmer_w_hits_masked_windows_wd_all %>% 
  filter(locus %in% endo_loci) %>%
  left_join(df_locus_order) %>% 
  mutate(locus = fct_reorder(locus, locus_order))

hits_TFBS_devCRE_all3 <- hits_TFBS_devCRE_all2 %>% 
  filter(locus %in% endo_loci) %>%
  left_join(df_locus_order) %>% 
  mutate(locus = fct_reorder(locus, locus_order))

df_quantiles2 <- df_quantiles%>% 
  filter(locus %in% endo_loci) %>%
  left_join(df_locus_order) %>% 
  mutate(locus = fct_reorder(locus, locus_order))


gata_hits <- hits_TFBS_devCRE_all3 %>% filter(!is.na(is_promoter) & category2=="specific, active") %>% pull(Gata4)
gata_all <- df_kmer_w_hits_masked_windows_wd_all2 %>% filter(!is.na(Gata4)) %>% pull(Gata4)

mean(gata_hits)/mean(gata_all)


plt_locis_TFBS <- ggplot(df_kmer_w_hits_masked_windows_wd_all2)+
  stat_bin2d(aes(x=Gata4, y=Sox17,colour=..count..),breaks=seq(-0.5,23))+guides(color="none")+ #max(df_kmer_w_hits_masked_windows_wd_all2$Gata4,na.rm=TRUE)
  geom_vline(data=df_quantiles2,aes(xintercept=q95_Gata4),color="coral1",linetype="dashed",size=0.3)+
  geom_vline(data=df_quantiles2,aes(xintercept=q99_Gata4),color="red",linetype="dashed",size=0.3)+
  new_scale_color()+
  new_scale_fill()+
  geom_point(data=hits_TFBS_devCRE_all3 %>% filter(!is.na(is_promoter)),aes(x=Gata4,y=Sox17,color=category2),size=1.5)+
  scale_color_manual(values=col_dict2)+
  facet_wrap(~locus,ncol=8)+
  coord_fixed()+
  labs(x="Number Gata4 TFBS [500 bp window]", y="Number Sox17 TFBS [500 bp window]") +
  theme_cowplot()+
  theme(panel.grid.major=element_line(size=0.1,color="grey"),
        strip.text=element_text(size=10),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        strip.background=element_blank())
plt_locis_TFBS
