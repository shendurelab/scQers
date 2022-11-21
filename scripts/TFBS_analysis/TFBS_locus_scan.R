
# lib dependencies
library(seqinr)
library(tidyverse)
library(Biostrings)
library(TFBSTools)
#library(catecolors)
library(pracma)
library(patchwork)
library(igraph)
library(universalmotif)

#library(ggnewscale)


args <- commandArgs(trailingOnly=TRUE)
gene_oi <- args[1]


dir_cluster <- "locus_level_scan_TFBS"
dir_cluster_dependencies <- "dependencies"

# scripts required for scanning 8-mer to fix orientation.
source(sprintf("%s/PWM_scan_seq_oi_motif_oi_v2_20220206.R",dir_cluster_dependencies))
source(sprintf("%s/oneHot_seq_20220124.R",dir_cluster_dependencies))
source(sprintf("%s/plot_sequence_with_motifs_v2_20220710.R",dir_cluster_dependencies))
source(sprintf("%s/PWM_scan_util_20220605.R",dir_cluster_dependencies))



# load and pre-process TFBS affinity table (empirically derived)
file_affinity_tables <- sprintf("%s/TFBS_relative_affinity_tables_endo_TFs_Uniprobe.txt",dir_cluster_dependencies)
df_affinity <- read.table(file_affinity_tables,header=TRUE)
df_affinity_long <- df_affinity %>% pivot_longer(cols=c(Gata4,Foxa2,Sox17))
df_affinity_long_for <- df_affinity_long %>% select(kmer_seq=for_8mer,TF_name=name,for_value=value)
df_affinity_long_rev <- df_affinity_long %>% select(kmer_seq=RC_8mer,TF_name=name,rev_value=value)


find_localMaxima <- function(x) {
  
  # dealing with NAs
  x[is.na(x)] <- 0
  
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]){ 
    y <- y[-1]
  }
  y
}



list_core_PWMs<- sprintf("%s/PFMs_TF_endo_core.txt",dir_cluster_dependencies)
pfms_combined <- readJASPARMatrix(list_core_PWMs, matrixClass=c("PFM"))
motif_names <- names(pfms_combined)
PWMs_list <- c()
for (motif in motif_names){
  pfm_oi <- pfms_combined[[motif]]
  PWMs_list <- c(PWMs_list,toPWM(pfm_oi))
}
PWMs <- do.call(PWMatrixList,PWMs_list)
PWMs_ids <- ID(PWMs)
names(PWMs) <- ID(PWMs)

get_TFBS_orientation <- function(df_hits,PWMs){
  
  TFBS_orientation <- c()
  
  for (idx in seq(dim(df_hits)[1])){
    df_PWM_scan <- PWM_scan_seq_oi_motif_oi_v2_20220206(df_hits$kmer_seq[idx] %>% DNAStringSet(),
                                                        PWMs[[df_hits$TF_name[idx]]])
    # get max value's orientation
    orientation_max <- df_PWM_scan %>% filter(score==max(score,na.rm=TRUE)) %>% pull(strand) %>% unique()
    if (length(orientation_max)>1){
      TFBS_orientation[idx] <- "0"
    } else {
      TFBS_orientation[idx] <- orientation_max
    }
  }
  
  return(TFBS_orientation)
}
get_collapsed_TFBS_hits_v2 <- function(df_kmer_w_affinity,
                                       affinity_threshold,
                                       dist_thresh){
  
  list_TFs <- unique(df_kmer_w_affinity %>% pull(TF_name))
  
  df_TFBS_hits <- data.frame()
  
  for (TF_oi in list_TFs){
    
    print(TF_oi)
    
    # restrict to TF of interest
    df_TF_oi <- df_kmer_w_affinity %>% filter(TF_name==TF_oi)
    
    # identify local maxima
    aff_trace <-df_TF_oi %>% pull(value)
    
    local_max <- find_localMaxima(aff_trace)
    # local_max <- which(diff(sign(diff(aff_trace)))==-2)+1
    
    # called all motifs positions above threshold
    pos_above_thresh <- df_TF_oi$start_pos[aff_trace>affinity_threshold[TF_oi]]
    
    
    # local maxima
    motif_oi_hits <- intersect(pos_above_thresh,local_max)
    df_TF_oi_hits <- df_TF_oi %>% filter(start_pos %in% motif_oi_hits)
    
    if (!isEmpty(df_TF_oi_hits$start_pos)){
      
      # get orientation of motifs
      TFBS_oi_orientation <- get_TFBS_orientation(df_TF_oi_hits,PWMs)
      
      df_TF_oi_hits2 <- df_TF_oi_hits %>% transform(TFBS_orientation=TFBS_oi_orientation)
      
      # loop through the different orientation
      for (ori in unique(df_TF_oi_hits2$TFBS_orientation)){
        
        df_TF_oi_hits2_ori <- df_TF_oi_hits2 %>% filter(TFBS_orientation==ori)
        
        # connected component of positions within a distance threshold
        dist_mat <- dist(df_TF_oi_hits2_ori$start_pos,diag = TRUE, upper = TRUE) %>% as.matrix
        conn_mat <- dist_mat<=dist_thresh
        g  <- graph.adjacency(conn_mat)
        clu <- components(g)
        
        # collapse positions above a threshold to the max position
        df_TF_oi_collapsed <- data.frame()
        for (mem in unique(clu$membership)){
          
          df_TF_oi_hits_cluster <- df_TF_oi_hits2_ori %>% filter(clu$membership==mem)
          df_TF_oi_hits_cluster_max <- df_TF_oi_hits_cluster %>% filter(value==max(value))
          
          df_TF_oi_collapsed <- rbind(df_TF_oi_collapsed,df_TF_oi_hits_cluster_max)
        }
        df_TFBS_hits <- rbind(df_TFBS_hits,df_TF_oi_collapsed)
      }
    }
  }
  
  return(df_TFBS_hits)
}
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}



# 200 kb sequences around TSS of considered loci
all_seqs <- readRDS(sprintf("%s/endo_loci_pm_100kb_20220915.RDS",dir_cluster_dependencies))

affinity_threshold <- c(0.3,0.125,0.3) # different thresholds selected from statistical properties of 8-mer relative affinities.
names(affinity_threshold) <- c("Gata4","Foxa2","Sox17")
dist_thresh <- 3
kmer_length <- 8
dy_ori <- 1


df_affinity_thershold <- data.frame(affinity_threshold,TF_name=names(affinity_threshold))

df_TFBS_hits_all_species <- data.frame()



seq_oi <- all_seqs[[gene_oi]]
length_seq <- str_length(seq_oi)


name_oi <- names(seq_oi)
print(name_oi)


# seq_char <- as.character(seq_oi)
seq_char_vec <- str_split(as.character(seq_oi),"") %>% unlist

df_kmer <- data.frame(matrix(nrow=length(seq(1,str_length(seq_oi)-kmer_length+1)),ncol=3))
colnames(df_kmer) <- c("start_pos","end_pos","kmer_seq")

for (i in seq(1,str_length(seq_oi)-kmer_length+1)){
  if ((i %% 5000)==0){
    print(i)
  }
  df_kmer[i,1] <- i
  df_kmer[i,2] <- i+kmer_length-1
  df_kmer[i,3] <- paste0(seq_char_vec[i:(i+kmer_length-1)],collapse="")
}

print("done with k-mer data frame construction")


# updated to deal with NAs
df_kmer_w_affinity <- data.frame()

for (TF_oi in unique(df_affinity_long_rev$TF_name)){
  
  print(TF_oi)
  
  df_kmer_w_affinity_oi <- df_kmer %>% left_join(df_affinity_long_for %>% filter(TF_name==TF_oi) )
  df_kmer_w_affinity2_oi <- df_kmer_w_affinity_oi %>% left_join(df_affinity_long_rev %>% filter(TF_name==TF_oi),by = "kmer_seq")
  
  TF_name <- df_kmer_w_affinity2_oi$TF_name.x
  TF_name[is.na(TF_name)] <- df_kmer_w_affinity2_oi$TF_name.y[is.na(TF_name)]
  
  value <- df_kmer_w_affinity2_oi$for_value
  value[is.na(value)] <- df_kmer_w_affinity2_oi$rev_value[is.na(value)]
  
  df_kmer_w_affinity3_oi <- df_kmer_w_affinity2_oi %>% 
    transform(TF_name=TF_name, value=value) %>%
    select(-TF_name.x,-TF_name.y,-for_value,-rev_value)
  
  df_kmer_w_affinity4_oi <- df_kmer_w_affinity3_oi
  df_kmer_w_affinity4_oi[is.na(df_kmer_w_affinity4_oi$TF_name),"TF_name"] <- TF_oi
  
  df_kmer_w_affinity <- rbind(df_kmer_w_affinity,df_kmer_w_affinity4_oi)
  
}

# as of 07/05/2022, not sure why this is here
df_kmer_w_affinity4 <- df_kmer_w_affinity %>% unique


# get local maxima for each trace
df_TFBS_hits <- get_collapsed_TFBS_hits_v2(df_kmer_w_affinity4,affinity_threshold,dist_thresh)


if (dim(df_TFBS_hits)[1]>0){
  
  df_TFBS_hits$TFBS_hit <- TRUE
  
  df_kmer_w_affinity5 <- df_kmer_w_affinity4 %>% left_join(df_TFBS_hits)
  df_kmer_w_affinity5$TFBS_hit[is.na(df_kmer_w_affinity5$TFBS_hit)] <- FALSE
  
  
  TF_ori_pos <- df_TFBS_hits$TFBS_orientation
  TF_ori_pos[TF_ori_pos=="+"] <- dy_ori
  TF_ori_pos[TF_ori_pos=="-"] <- -dy_ori
  TF_ori_pos[TF_ori_pos=="0"] <- 0
  TF_ori_pos <- as.numeric(TF_ori_pos)
  #   
  df_TFBS_hits_w_pos <- df_TFBS_hits %>% transform(TF_ori_pos=TF_ori_pos)
  df_TFBS_hits_w_pos2 <- df_TFBS_hits_w_pos %>% transform(sequence_name=name_oi,
                                                          gene=gene_oi)
  #   
  df_TFBS_hits_all_species <- rbind(df_TFBS_hits_all_species,df_TFBS_hits_w_pos2)
  
}

# }

write.table(df_TFBS_hits_all_species,sprintf("TFBS_hits_%s_100kb_pm_TSS.txt",gene_oi),
            sep="\t",row.names=FALSE,quote=FALSE)

