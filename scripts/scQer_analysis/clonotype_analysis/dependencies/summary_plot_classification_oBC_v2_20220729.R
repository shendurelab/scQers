
summary_plot_classification_oBC_v2_20220729 <- function(df_cBC_metadata_oBC,oBC_clonotypes,
                                                     n_oBC_thresh, p_to_top_thresh, recovered_oBC_frac_thresh,
                                                     path_out,run_name,sample_name){
  
  
  date_str <- str_replace_all(Sys.Date(),"-","") 

  
  plt2 <- ggplot(df_cBC_metadata_oBC)+geom_hex(aes(y=frac_oBC_assigned_to_top_clone,x=frac_oBC_recovered_top_clone))+
    annotate("rect", xmin = recovered_oBC_frac_thresh, xmax = 1.03, 
             ymin = p_to_top_thresh, ymax = 1.03,alpha = .1,fill = "red")+
    annotate("rect", xmin = 0, xmax = recovered_oBC_frac_thresh, 
             ymin = 0, ymax = 1.03,alpha = .1,fill = "cyan")+
    annotate("rect", xmin = recovered_oBC_frac_thresh, xmax = 1.03, 
             ymin = 0, ymax = p_to_top_thresh,alpha = .1,fill = "goldenrod1")+
    labs(x="Fraction of oBC detected from top clonotype ",y="Fraction of assigned oBC to top clonotype")+coord_fixed()
  
  

  
  # classification 
  
  cBC_likely_singlet <- df_cBC_metadata_oBC %>% filter(number_assigned_oBC>=n_oBC_thresh & 
                                                         frac_oBC_assigned_to_top_clone>p_to_top_thresh & 
                                                         frac_oBC_recovered_top_clone> recovered_oBC_frac_thresh) %>% pull(cBC)
  
  cBC_likely_missed_clonotype <- df_cBC_metadata_oBC %>% filter(frac_oBC_recovered_top_clone<=recovered_oBC_frac_thresh) %>% pull(cBC)
  
  cBC_likely_doublet <- df_cBC_metadata_oBC %>% filter(frac_oBC_recovered_top_clone>recovered_oBC_frac_thresh & 
                                                         frac_oBC_assigned_to_top_clone<=p_to_top_thresh) %>% pull(cBC)
  
  df_cBC_metadata_oBC2 <- df_cBC_metadata_oBC
  df_cBC_metadata_oBC2$category <- NA
  df_cBC_metadata_oBC2$category[df_cBC_metadata_oBC2$cBC %in% cBC_likely_singlet] <- "singlet"
  df_cBC_metadata_oBC2$category[df_cBC_metadata_oBC2$cBC %in% cBC_likely_doublet] <- "doublet"
  df_cBC_metadata_oBC2$category[df_cBC_metadata_oBC2$cBC %in% cBC_likely_missed_clonotype] <- "missed_clonotype"
  
  
  plt1 <- ggplot(df_cBC_metadata_oBC2)+geom_hex(aes(x=number_assigned_oBC,y=frac_oBC_assigned_to_top_clone),bins=50)+scale_x_log10()+
    geom_segment(aes(x=n_oBC_thresh,xend=n_oBC_thresh,y=p_to_top_thresh,yend=1.05),color="red",linetype="dashed",size=0.3)+
    geom_segment(aes(x=n_oBC_thresh,xend=150,y=p_to_top_thresh,yend=p_to_top_thresh),color="red",linetype="dashed",size=0.3)+
    labs(x="Number assigned oBC per cellBC (>30 UMI)",y="Fraction of assigned oBC to top clonotype")
  
  plt1_v2 <- ggplot(df_cBC_metadata_oBC2)+geom_hex(aes(x=number_assigned_oBC,y=frac_oBC_assigned_to_top_clone),bins=50)+scale_x_log10()+
    geom_segment(aes(x=n_oBC_thresh,xend=n_oBC_thresh,y=p_to_top_thresh,yend=1.05),color="red",linetype="dashed",size=0.3)+
    geom_segment(aes(x=n_oBC_thresh,xend=150,y=p_to_top_thresh,yend=p_to_top_thresh),color="red",linetype="dashed",size=0.3)+facet_wrap(~category)+
    labs(x="Number assigned oBC per cellBC (>30 UMI)",y="Fraction of assigned oBC to top clonotype")
  
  plt_design <- "
  #CD
  ###
  #A#
  BBB"
  
  singlet_cBC_to_clonotype <- df_cBC_metadata_oBC2 %>% filter(category=="singlet")
  
  df_n_assigned_per_clonotype <- singlet_cBC_to_clonotype %>% group_by(top_oBC_clone) %>% summarize(n_assigned_cellBC=length(cBC))
  # df_n_assigned_per_clonotype <- df_n_assigned_per_clonotype %>% rename(c("top_oBC_clone"="clonotype_id")) %>% data.frame
  df_n_assigned_per_clonotype <- df_n_assigned_per_clonotype %>% rename(c("clonotype_id"="top_oBC_clone")) %>% data.frame
  
  
  num_oBC_per_clone <- sort(sapply(oBC_clonotypes,length),decreasing=TRUE)
  quantile(num_oBC_per_clone,seq(0,1,0.1))
  
  
  # list of oBC per clonotype
  df_oBC_per_clonotype <- data.frame()
  
  for (clone_oi in names(oBC_clonotypes)){
    list_cBC2 <- df_cBC_metadata_oBC2 %>% filter(category=="singlet" & top_oBC_clone==clone_oi) %>% pull(cBC) %>% paste(collapse=",")
    # reps <- str_split(list_cBC,"_") %>% lapply("[[",1) %>% unlist
    # cBC <- str_split(list_cBC,"_") %>% lapply("[[",2) %>% unlist
    # list_cBC2 <- paste0(reps,"_",cBC) %>% paste(collapse=",")
    
    list_oBC <- oBC_clonotypes[[clone_oi]] %>% paste(collapse=",")
    df_oBC_per_clonotype <- rbind(df_oBC_per_clonotype,data.frame(clonotype_id=clone_oi,
                                                                  list_oBC,
                                                                  list_cBC2))
    
  }
  
  df_clonotype_metadata <- data.frame(clonotype_id = names(num_oBC_per_clone),
                                      MOI=num_oBC_per_clone)
  
  df_clonotype_metadata2<- df_clonotype_metadata %>% left_join(df_n_assigned_per_clonotype) %>% left_join(df_oBC_per_clonotype)
  
  plt3 <- ggplot(df_clonotype_metadata2) + geom_hex(aes(x=MOI,y=n_assigned_cellBC,color=..count..))+guides(color="none")+scale_x_log10()+scale_y_log10()
  
  table_name <- sprintf("%s/oBC_clonotypes_metadata_%s_%s_%s.txt",
                                    path_out,sample_name,run_name,date_str)
  write.table(df_clonotype_metadata2,table_name,sep="\t",quote=FALSE, row.names=FALSE)

  fig_name <- sprintf("%s/summary_fig_assignment_oBC_clonotypes_%s_%s_%s.pdf",
                      path_out,run_name,sample_name,date_str)
  pdf(fig_name,height=15,width=15)
  print(plt1+plt1_v2+plt2+plt3+plot_layout(design=plt_design))
  dev.off()
  
  
  
  df_clone_meta_info <- data.frame()
  
  for (clone_oi in names(oBC_clonotypes)){
    
  }
  
  
  return(df_cBC_metadata_oBC2)
  
  
}