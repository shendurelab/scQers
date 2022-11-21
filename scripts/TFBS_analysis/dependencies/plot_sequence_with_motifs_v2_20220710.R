
plot_sequence_with_motifs_v2_20220710 <- function(seq_region_oi,df_pos_hits,motif_color_map){
  
  
  df_pos_hits <- df_pos_hits %>% select(motif=TF_name,
                                    start=start_pos,
                                    end=end_pos,
                                    mostif_hit_seq=kmer_seq,
                                    strand=TFBS_orientation,
                                    score=value)

  
  # convert to long form, assumes all 8-mers
  df_pos_hits_long <- data.frame()
  
  for (hit_id in seq(dim(df_pos_hits)[1])){
    
    pos_range_oi <- seq(df_pos_hits$start[hit_id],df_pos_hits$end[hit_id])
    
    if (df_pos_hits$strand[hit_id] == "0"){
      df_pos_hits_long <- rbind(df_pos_hits_long,
                                data.frame(motif=df_pos_hits$motif[hit_id],
                                           positions=pos_range_oi,
                                           strand="+",
                                           score=df_pos_hits$score[hit_id]))
      df_pos_hits_long <- rbind(df_pos_hits_long,
                                data.frame(motif=df_pos_hits$motif[hit_id],
                                           positions=pos_range_oi,
                                           strand="-",
                                           score=df_pos_hits$score[hit_id]))
      
    } else {
      df_pos_hits_long <- rbind(df_pos_hits_long,
                                data.frame(motif=df_pos_hits$motif[hit_id],
                                           positions=pos_range_oi,
                                           strand=df_pos_hits$strand[hit_id],
                                           score=df_pos_hits$score[hit_id]))
    }
    
  }
  
  # # # # # # # # # # # # # # # # # # # # 
  # function to plot sequence w/ motif hits
  # # # # # # # # # # # # # # # # # # # # 
  colScale <- scale_colour_manual(name = "motif",values = motif_color_map)
  
  # fixed parameters for aesthetics. 
  dy_strand <- 0.02
  dy_break <- 0.1
  font_size <- 3
  line_break <- 150
  
  seq_region_oi <- seq_region_oi %>% DNAStringSet()
  length_seq <- width(seq_region_oi)
  df_text_seq_region <- data.frame(letters=c(strsplit(as.character(seq_region_oi),"")[[1]], strsplit(as.character(complement(seq_region_oi)),"")[[1]]),
                                   positions= c(seq(length_seq), seq(length_seq)), 
                                   y= c(rep(0,length_seq), rep(-dy_strand,length_seq)),
                                   strand=c(rep("+",length_seq), rep("-",length_seq)))
  
  
  # add the motif hit annotation: 
  x_break <- df_text_seq_region$positions %% line_break
  x_break[x_break==0] <- line_break
  y_break <- -dy_break*floor((df_text_seq_region$positions-1)/line_break)
  y_break <- y_break+df_text_seq_region$y
  df_text_seq_region$x_break <- x_break
  df_text_seq_region$y_break <- y_break
  
  
  if (dim(df_pos_hits)[1]>0){
    df_text_seq_region2 <- left_join(df_text_seq_region,df_pos_hits_long)
    df_text_seq_region2$font_type <- "plain"
    df_text_seq_region2$font_type[!is.na(df_text_seq_region2$motif)] <- "bold"
  } else {
    df_text_seq_region2$motif <- NA
    df_text_seq_region2$font_type <- "plain"
  }
  
  plt_seq_motif <- ggplot(df_text_seq_region2)+
    geom_text(aes(x=x_break,y=y_break,label=letters,color=motif,fontface=font_type),family="Courier",size=font_size)+theme_void()+colScale 
  
  
  
}
