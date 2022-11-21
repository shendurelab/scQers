library(Seurat)
library(dplyr)
library(patchwork)
#library(sctransform)
#library(ggplot2)
library(tidyverse)
library(scales)
#library(tictoc)
#library(irlba) 
#library(schex) # might need to install
#library(stringdist)
library(Biostrings)
library(stringr)
library(cowplot)


# load in the data
df_oBC <- read.table("table_oBC_counts_w_gex_all_reps_v2_all_umis_gex_valid_20220824.txt",
                     header=TRUE)

GEx_obj_final <-readRDS("/Users/jbl/Documents/UW/manuscripts/single_cell_expression_reporter/supp_tables/GEx_mEB_final_20221025.RDS")
valid_cells <- Cells(GEx_obj_final)
df_oBC <- df_oBC %>% filter(full_cellBC_pdT %in% valid_cells) %>% transform(biol_rep=substr(rep_id,1,str_length(rep_id)-1))

df_oBC2 <- df_oBC %>% filter(p25_valid_oBC & CRE_class=="devCRE") #%>% left_join(df_mBC_pdT) %>% select(-c(capture_modality,mBC_in_full_list))



# join in clone information: see generate_final_clonotype_lists.R
df_clone_assignment <- read.table("metadata_cell_final_clonotype_assignments_all_reps_mEB_20221027.txt",
                                  header=TRUE)
clone_oBC_list <- readRDS("list_final_oBC_clonotypes_all_reps_mEB_20221027.RDS")

df_singlets <- df_clone_assignment %>% filter(category=="singlet") %>% select(full_cellBC_pdT=cBC, clone_id=top_oBC_clone)
df_oBC_singlets <- df_oBC2 %>% left_join(df_singlets) %>% filter(!is.na(clone_id))



# build a data frame averaging over the single-cell data for clones: 
clone_ids <- names(clone_oBC_list)
count_thresh <- seq(0,15)
n_cell_per_clone_thresh <- 2
df_oBC_dropout <- c()

for (clone_id_oi in clone_ids){
  
  
  clone_oBC_list_oi <- clone_oBC_list[[clone_id_oi]]
  
  sc_rep_mBC_oBC_clone_oi <- df_oBC_singlets %>% filter(clone_id==clone_id_oi) 
  
  cBC_clone <- unique(sc_rep_mBC_oBC_clone_oi$full_cellBC_pdT)
  n_cells_per_clone <- length(cBC_clone)
  
  if (length(cBC_clone)>n_cell_per_clone_thresh){
    
    print(clone_id_oi)
    
    # not restricting attention to oBC uniquely paired to mBC for this analysis
    sc_rep_mBC_oBC_clone_oi2 <- sc_rep_mBC_oBC_clone_oi %>% select(full_cellBC_pdT,oBC,reads_oBC,UMIs_oBC) %>% unique()
    
    n_oBC_in_clone <- length(unique(sc_rep_mBC_oBC_clone_oi2$oBC))
    
    for (idx in seq(length(count_thresh))){
      
      df_dropout_oi <-  sc_rep_mBC_oBC_clone_oi2 %>% group_by(full_cellBC_pdT) %>% 
        summarize( TP= sum( (oBC %in% clone_oBC_list_oi) & (UMIs_oBC>count_thresh[idx])),
                   FP= sum( !(oBC %in% clone_oBC_list_oi) & (UMIs_oBC>count_thresh[idx])),
                   FN= length(clone_oBC_list_oi)-sum( (oBC %in% clone_oBC_list_oi) & (UMIs_oBC>count_thresh[idx]))) %>%
                   # TN= n_oBC_in_clone- 
                   #   sum( !(oBC %in% clone_oBC_list_oi) & (UMIs_oBC>count_thresh[idx]))-
                   #   length(clone_oBC_list_oi)) %>%
        data.frame() %>% transform(n_cell_in_clone=n_cells_per_clone)
      
      df_dropout_oi$n_oBC_clone <- length(clone_oBC_list_oi)
      df_dropout_oi$oBC_UMI_thresh <- count_thresh[idx]
      df_dropout_oi$clone_id <- clone_id_oi
      df_dropout_oi$cell_line <- sc_rep_mBC_oBC_clone_oi$gex_cluster_id[1]
      df_dropout_oi$rep_id <- sc_rep_mBC_oBC_clone_oi$rep_id[1]
      
      df_oBC_dropout <- rbind(df_oBC_dropout,df_dropout_oi)
      
    }
  }
}





df_oBC_dropout2 <- df_oBC_dropout %>% transform(biol_rep=str_sub(rep_id,1,str_length(rep_id)-1))

thresh_n_cells <- 3
df_oBC_dropout2 <- df_oBC_dropout %>% filter(n_cell_in_clone>=thresh_n_cells) %>% transform(biol_rep=str_sub(rep_id,1,str_length(rep_id)-1))


# aggregate quantificaiton 
df_summed_dropout <- df_oBC_dropout2 %>% group_by(oBC_UMI_thresh,biol_rep) %>% 
  summarize(total_num_FP=sum(FP), 
            total_num_FN=sum(FN),
            total_num_TP=sum(TP),
            true_positive_rate=sum(TP)/(sum(TP)+sum(FN)),
            false_discovery_rate=sum(FP)/(sum(FP)+sum(TP)),
            false_negative_rate=sum(FN)/(sum(TP)+sum(FN)),
            precision=sum(TP)/(sum(TP)+sum(FP)),
            recall=sum(TP)/(sum(TP)+sum(FN)),
            total_num_oBCs=sum(n_oBC_clone)) 



# summary plots
plt_PR <- ggplot(df_summed_dropout)+
  geom_step(aes(x=precision,y=recall,color=biol_rep))+
  coord_fixed(ylim=c(0.94,1),xlim=c(0.96,1))+
 
  guides(color="none")+theme_cowplot()+
  theme(panel.grid.major=element_line(size=0.1,color="grey"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))+
  labs(x="Precision",
       y="Recall")
plt_PR

plt3 <- ggplot(df_summed_dropout)+
  geom_step(aes(x=oBC_UMI_thresh+1,y=false_discovery_rate,color=biol_rep))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(1E-3,1))+
  xlim(c(0,13))+
  guides(color="none")+annotation_logticks(sides="l")+theme_cowplot()+
  theme(panel.grid.major=element_line(size=0.1,color="grey"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))+
  labs(x="oBC UMI threshold",y="False discovery rate\n[FP/(TP+FP)]")
plt3_2 <- ggplot(df_summed_dropout)+
  geom_step(aes(x=oBC_UMI_thresh+1,y=false_negative_rate,col=biol_rep))+
  xlim(c(0,13))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(1E-3,0.1))+
  annotation_logticks(sides="l")+theme_cowplot()+
  theme(panel.grid.major=element_line(size=0.1,color="grey"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))+
  labs(x="oBC UMI threshold",y="False negative rate\n[FN/(TP+FN)]")




