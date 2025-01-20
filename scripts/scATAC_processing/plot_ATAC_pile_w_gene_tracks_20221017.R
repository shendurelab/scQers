
plot_ATAC_pile_w_gene_tracks_20221017 <- function(df_peak_set,
                                                  proj_QC_no_doub,
                                                  extension_TSS,
                                                  locus_oi,lineage_all,
                                                  df_DE_all,
                                                  df_reporter_results2,
                                                  bool_all_clusters,
                                                  annotation_for_tracks){
  
  library(patchwork)
  
  geneAnnotation <-  getGeneAnnotation(proj_QC_no_doub)
  regions <- geneAnnotation$genes
  
  lineage_oi <- lineage_all %>% filter(DE_gene==locus_oi) %>% pull(DE_lineage)
  print(locus_oi)
  
  
  region_oi <- regions[which(tolower(mcols(regions)$symbol) %in% tolower(locus_oi))]
  region_oi <- region_oi[order(match(tolower(mcols(region_oi)$symbol), tolower(locus_oi)))]
  region_oi_TSS <- resize(region_oi, 1, "start")
  strand(region_oi_TSS) <- "*"
  peak_oi_range_ext_long2 <- extendGR(region_oi_TSS, upstream = extension_TSS, downstream = extension_TSS)
  
  
  
  chr_pos_str <- paste0(peak_oi_range_ext_long2@seqnames[1] %>% as.character,
                        ":",peak_oi_range_ext_long2@ranges@start,
                        "-",peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width, " (mm10)")
  
  max_cells <- 1E5
  size_tile <- 150
  df <- ArchR:::.groupRegionSumArrows(ArchRProj = proj_QC_no_doub, groupBy = "Clusters_coarse", 
                                      normMethod = "ReadsInTSS", useGroups = NULL, minCells = 400, 
                                      region = peak_oi_range_ext_long2, tileSize = size_tile, threads = 4, 
                                      verbose = TRUE, maxCells=max_cells)
  df <- df %>% filter(group %in% c("endo","meso","ecto","pluri"))
  
  
  
  
  dev_CREs_loci <- df_DE_all %>% filter(gene_class=="devCRE_locus") %>% pull(DE_gene) %>% unique()
  df_DE_all_locus_oi <- df_DE_all %>% filter(DE_gene==locus_oi & DE_test_in_germ_layer==lineage_oi )
  germ_layer_oi <- df_DE_all_locus_oi$DE_lineage[1]
  targeted_peaks <- df_DE_all_locus_oi %>% pull(CRE_id) %>% unique
  
  df_peak_set_oi <- df_peak_set %>% filter(unique_peak_id %in% targeted_peaks)
  
  
  
  df_peak_set_oi2 <- df_peak_set_oi %>% 
    left_join(df_reporter_results2 %>% 
                dplyr:::select(unique_peak_id,
                               all_rep_act_hit,
                               all_rep_spec_hit,
                               is_promoter)) %>% 
    transform(group="ecto")
  
  buffer_peak_highlight <- 200
  
  # col_groups <- c("ecto"="#66C2A5",
  #                 "endo"="#FC8D62",
  #                 "meso"="#8DA0CB",
  #                 "pluri"="404040")
  
  col_groups <- c("ecto"="#66C2A5",
                  "endo"="#FC8D62",
                  "meso"="#8DA0CB",
                  "pluri"="#E78AC3")
  
  
  df_peak_set_oi3 <- df_peak_set_oi %>% 
    left_join(df_reporter_results2 %>% 
                dplyr:::select(unique_peak_id,
                               all_rep_act_hit,
                               all_rep_spec_hit,
                               is_promoter)) %>% 
    transform(group=lineage_oi)
  
  buffer_peak_highlight <- 250
  
  margin_dx <- 0
  
  
  if (sum(df_peak_set_oi2$all_rep_act_hit,na.rm=TRUE)>0){
    plt_ATAC2 <- ggplot() + 
      geom_rect(data=df_peak_set_oi,
                aes(xmin=start-buffer_peak_highlight,xmax=end+buffer_peak_highlight,ymin=0,ymax=max(df$y)),
                fill="#E6E6E6")+
      geom_bar(data=df,aes(x=x,y=y,fill=group,color=group),size=0.46,stat = "identity")+
      geom_point(data=df_peak_set_oi2 %>% filter(all_rep_act_hit & is_promoter),inherit.aes=FALSE,
                 aes(x=start+250,y=0.9*max(df$y)),size=2,shape=25,fill="goldenrod",color="black",stroke=0.4)+
      geom_point(data=df_peak_set_oi2 %>% filter(all_rep_act_hit & !is_promoter & !all_rep_spec_hit),inherit.aes=FALSE,
                 aes(x=start+250),y=0.9*max(df$y),size=2,shape=25,fill="black",color="black",stroke=0.4)+
      geom_point(data=df_peak_set_oi2 %>% filter(all_rep_spec_hit),inherit.aes=FALSE,
                 aes(x=start+250),y=0.9*max(df$y),size=2,shape=25,fill="red",color="black",stroke=0.4)+
      scale_color_manual(values=col_groups)+
      scale_fill_manual(values=col_groups)+
      scale_x_continuous(limits= c(peak_oi_range_ext_long2@ranges@start,peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width),
                         expand=expansion(c(0,0)))+
      facet_wrap(facets = ~group,
                 ncol = 1)+
      guides(fill="none",color="none")+theme_void()+
      theme(plot.caption =element_text(size=5,hjust = 1),
            strip.text.x = element_blank(),
            plot.margin = unit(x=c(margin_dx,-margin_dx,0.1,-margin_dx),units="inches"))+
            # panel.border=element_rect(color="grey",size=0.2))+
      labs(caption=sprintf("%s, %s",locus_oi,chr_pos_str))
    
    
    
    
    # plt1+theme_void()
    
    
    plt_ATAC_two_lineages <- ggplot() + 
      geom_rect(data=df_peak_set_oi, inherit.aes=FALSE,
                aes(xmin=start-buffer_peak_highlight,xmax=end+buffer_peak_highlight,ymin=0,ymax=max(df$y)),
                fill="#E6E6E6")+
      geom_bar(data=df %>% filter(group %in% c("pluri",lineage_oi)),
               aes(x=x,y=y,fill=group,color=group),size=0.5,stat = "identity")+
      geom_point(data=df_peak_set_oi3 %>% filter(all_rep_act_hit & is_promoter),inherit.aes=FALSE,
                 aes(x=start+250,y=0.9*max(df$y)),size=2,shape=25,fill="goldenrod",color="black",stroke=0.4)+
      geom_point(data=df_peak_set_oi3 %>% filter(all_rep_act_hit & !is_promoter),inherit.aes=FALSE,
                 aes(x=start+250),y=0.9*max(df$y),size=2,shape=25,fill="black",color="black",stroke=0.4)+
      geom_point(data=df_peak_set_oi3 %>% filter(all_rep_spec_hit),inherit.aes=FALSE,
                 aes(x=start+250),y=0.9*max(df$y),size=2,shape=25,fill="red",color="black",stroke=0.4)+
      scale_color_manual(values=col_groups)+
      scale_x_continuous(limits= c(peak_oi_range_ext_long2@ranges@start,peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width),
                         expand=expansion(c(0,0)))+
      facet_wrap(facets = ~group, ncol = 1)+
      theme_void()+guides(fill="none",color="none")+
      theme(plot.caption =element_text(size=5,hjust = 1),
            strip.text.x = element_blank(),
            plot.margin = unit(x=c(margin_dx,-margin_dx,0.1,-margin_dx),units="inches"))+
      labs(caption=chr_pos_str)
    
  } else {
    plt_ATAC2 <- ggplot() + 
      geom_rect(data=df_peak_set_oi, inherit.aes=FALSE,
                aes(xmin=start-buffer_peak_highlight,xmax=end+buffer_peak_highlight,ymin=0,ymax=max(df$y)),
                fill="#E6E6E6")+
      geom_bar(data=df, aes(x=x,y=y,fill=group,color=group),size=0.5,stat = "identity")+
      scale_color_manual(values=col_groups)+
      facet_wrap(facets = ~group,
                 ncol = 1)+
      scale_x_continuous(limits= c(peak_oi_range_ext_long2@ranges@start,peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width),
                         expand=expansion(c(0,0)))+
      theme_void()+guides(fill="none",color="none")+
      theme(plot.caption =element_text(size=5,hjust = 1),
            strip.text.x = element_blank(),
            plot.margin = unit(x=c(margin_dx,-margin_dx,0.1,-margin_dx),units="inches"))+
      labs(caption=sprintf("%s, %s",locus_oi,chr_pos_str))

    plt_ATAC_two_lineages <- ggplot() + 
      geom_rect(data=df_peak_set_oi, inherit.aes=FALSE,
                aes(xmin=start-buffer_peak_highlight,xmax=end+buffer_peak_highlight,ymin=0,ymax=max(df$y)),
                fill="#E6E6E6")+
      geom_bar(data=df %>% filter(group %in% c("pluri",lineage_oi)),
               aes(x=x,y=y,fill=group,color=group),size=0.5,stat = "identity")+
      scale_color_manual(values=col_groups)+
      scale_x_continuous(limits= c(peak_oi_range_ext_long2@ranges@start,peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width),
                         expand=expansion(c(0,0)))+
      facet_wrap(facets = ~group, ncol = 1)+
      theme_void()+guides(fill="none",color="none")+
      theme(plot.caption =element_text(size=5,hjust = 1),
            strip.text.x = element_blank(),
            plot.margin = unit(x=c(margin_dx,-margin_dx,0.1,-margin_dx),units="inches"))+
      labs(caption=sprintf("%s, %s",locus_oi,chr_pos_str))
  }
  
  
  # peak_oi_range_ext_long_for_anno <- extendGR(region_oi_TSS, upstream = 2*extension_TSS, downstream = 2*extension_TSS)
  # 
  # 
  # plt_gene_track <- autoplot(annotation_for_tracks, 
  #                            which = peak_oi_range_ext_long_for_anno,
  #                            label=FALSE,
  #                            color = "grey",
  #                            fill = "grey",
  #                            coord = c("genome"),
  #                            size=0.3)+
  #   coord_cartesian(xlim= c(peak_oi_range_ext_long2@ranges@start,peak_oi_range_ext_long2@ranges@start+peak_oi_range_ext_long2@ranges@width),
  #                   expand=FALSE)+theme_void()
  
  
  if (bool_all_clusters){
    # plt_design <- "A
    # A
    # A
    # A
    # A
    # A
    # B"
    # return(plt_ATAC2+plt_gene_track@ggplot+plot_layout(design=plt_design))
    
    
    return(plt_ATAC2)
    
  } else {
    # plt_design <- "A
    # A
    # A
    # B"
    # return(plt_ATAC_two_lineages+plt_gene_track@ggplot+plot_layout(design=plt_design))
    
    return(plt_ATAC_two_lineages)
  }
  
}
