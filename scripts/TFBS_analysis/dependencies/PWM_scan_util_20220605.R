PWM_scan_util_20220605 <- function(file_seqs,file_PWMs,score_threshold,plot_bool,print_table_bool){
  
  # inputs:
  # file_seqs: path to sequences in fasta
  # file_PWMs: path to PWMs (JASPAR format)
  # score_threshold: log2-odds score threshold for calling a match, see Pizzi et al, "Finding significant matches of position weight matrices in linear time" 2011.  as an example ref
  # plot_bool: boolean whether to generate a pdf of sequence with matched motifs
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # function to convert to one-hot encoding
  # # # # # # # # # # # # # # # # # # # # 
  oneHot_seq_20220124 <- function(seq_oi, n = c("A", "C", "G", "T")) {
    
    # Construct matrix
    seq_oi <- toupper(seq_oi)
    seq_len <- nchar(seq_oi)
    seq_split <- unlist(strsplit(x = seq_oi, split = ""))
    seq_mat <- matrix(data = rep(0, seq_len * length(n)), nrow = 4)
    rownames(seq_mat) <- n
    colnames(seq_mat) <- seq_split
    
    # Encode
    for (i in n) {
      seq_mat[rownames(seq_mat) == i, colnames(seq_mat) == i] <- 1
    }
    return(seq_mat)
  }
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # function to scan PWMs across sequence
  # # # # # # # # # # # # # # # # # # # # 
  PWM_scan_seq_oi_motif_oi_v2_20220206 <- function(seq_oi,PWM_oi){
    
    # get reverse complement and convert to one-hot encoding
    seq_oi_RC <- reverseComplement(seq_oi)
    seq_size <- width(seq_oi)
    oh_seq <- Reshape(oneHot_seq_20220124(seq_oi),4,seq_size)
    oh_seq_RC <- Reshape(oneHot_seq_20220124(seq_oi_RC),4,seq_size)
    
    name_motif <- name(PWM_oi)
    PWM_oi <- Matrix(PWM_oi)
    region_size <- length(oh_seq[1,])
    
    # convolution on both strands
    convoA <- stats:::filter(oh_seq[1,],rev(PWM_oi[1,]), method = "convolution", sides=1) %>% as.numeric()
    convoC <- stats:::filter(oh_seq[2,],rev(PWM_oi[2,]), method = "convolution", sides=1) %>% as.numeric()
    convoG <- stats:::filter(oh_seq[3,],rev(PWM_oi[3,]), method = "convolution", sides=1) %>% as.numeric()
    convoT <- stats:::filter(oh_seq[4,],rev(PWM_oi[4,]), method = "convolution", sides=1) %>% as.numeric()
    convo <- convoA+convoC+convoG+convoT
    
    convoA_RC <- stats:::filter(oh_seq_RC[1,],rev(PWM_oi[1,]), method = "convolution", sides=1) %>% as.numeric()
    convoC_RC <- stats:::filter(oh_seq_RC[2,],rev(PWM_oi[2,]), method = "convolution", sides=1) %>% as.numeric()
    convoG_RC <- stats:::filter(oh_seq_RC[3,],rev(PWM_oi[3,]), method = "convolution", sides=1) %>% as.numeric()
    convoT_RC <- stats:::filter(oh_seq_RC[4,],rev(PWM_oi[4,]), method = "convolution", sides=1) %>% as.numeric()
    convo_reverse <- rev(convoA_RC+convoC_RC+convoG_RC+convoT_RC)
    
    # collect to a data frame and bind
    df_PWM_oi_for <- data.frame(relative_positions=1:length(convo), score=convo, strand="+", motif=name_motif)
    df_PWM_oi_rev <- data.frame(relative_positions=1:length(convo), score=convo_reverse, strand="-", motif=name_motif)
    df_PWM_oi_max <- data.frame(relative_positions=1:length(convo), score=pmax(convo,convo_reverse), strand="max", motif=name_motif)
    df_PWM_scores <- rbind(df_PWM_oi_for,df_PWM_oi_rev,df_PWM_oi_max)
    
    return(df_PWM_scores)
  }
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # function to get motif hit from PWM scan
  # # # # # # # # # # # # # # # # # # # # 
  get_position_motif_hit_df_no_offset_20220207 <- function(df_PWM_scores,name_seq_oi,thresh_oi,PWMs,sequence){
    
    motif_hits <- df_PWM_scores %>% filter(seq==name_seq_oi & score>thresh_oi & strand!="max")

    df_pos_hits_short <- data.frame()
    df_pos_hits_long <- data.frame()
    
    if (dim(motif_hits)[1]>0){
      for (ind in seq(dim(motif_hits)[1])){
        
        motif_length <- dim(PWMs[[motif_hits$motif[ind]]])[2]
        
        if (motif_hits$strand[ind]=="+"){
          start_oi <- -motif_length+1 #df_motif_range %>% filter(motif==motif_hits$motif[ind] & strand==motif_hits$strand[ind]) %>% pull(start)
          end_oi <- 0 
        } else {
          start_oi <- 0 
          end_oi <- motif_length-1
        }
        
        pos_oi <- motif_hits$relative_positions[ind]
        
        #browser()
        
        if (!isempty(start_oi)){
          pos_range_oi <- pos_oi+start_oi:end_oi
          
          # tryCatch(rbind(df_pos_hits,data.frame(motif=motif_hits$motif[ind],
          #                                       positions=pos_range_oi,
          #                                       strand=motif_hits$strand[ind])),
          #          finally=browser())
          
          df_pos_hits_short <- rbind(df_pos_hits_short,data.frame(motif=motif_hits$motif[ind],
                                                      start=pos_oi+start_oi,
                                                      end=pos_oi+end_oi,
                                                      motif_hit_seq=substr(sequence,min(pos_range_oi),max(pos_range_oi)),
                                                      motif_hit_seq_RC=substr(sequence,min(pos_range_oi),max(pos_range_oi)) %>% DNAStringSet() %>% reverseComplement,
                                                      strand=motif_hits$strand[ind],
                                                      score=motif_hits$score[ind]))
          
          df_pos_hits_long <- rbind(df_pos_hits_long,data.frame(motif=motif_hits$motif[ind],
                                                      positions=pos_range_oi,
                                                      strand=motif_hits$strand[ind],
                                                      score=motif_hits$score[ind]))
          
        }
      }
    }
    return(list(df_pos_hits_short,df_pos_hits_long))
  }
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # function to plot sequence w/ motif hits
  # # # # # # # # # # # # # # # # # # # # 
  plot_sequence_with_motifs_20220202 <- function(seq_region_oi,df_pos_hits,motif_color_map){
    
    colScale <- scale_colour_manual(name = "motif",values = motif_color_map)
    
    # fixed parameters for aesthetics. 
    dy_strand <- 0.02
    dy_break <- 0.1
    font_size <- 5
    line_break <- 200
    
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
      df_text_seq_region2 <- left_join(df_text_seq_region,df_pos_hits)
      df_text_seq_region2$font_type <- "plain"
      df_text_seq_region2$font_type[!is.na(df_text_seq_region2$motif)] <- "bold"
    } else {
      df_text_seq_region2$motif <- NA
      df_text_seq_region2$font_type <- "plain"
    }
    
    plt_seq_motif <- ggplot(df_text_seq_region2)+
      geom_text(aes(x=x_break,y=y_break,label=letters,color=motif,fontface=font_type),family="Courier",size=font_size)+theme_void()+colScale 
    return(plt_seq_motif)

  }
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # read in sequences (assumed fasta)
  # # # # # # # # # # # # # # # # # # # # 
  print("reading seqs")
  seqs <- read.fasta(file_seqs, as.string = TRUE)
  name_seqs <- names(seqs)
  seqs <- seqs %>% as.character()
  names(seqs) <- name_seqs
  seqs <- seqs %>% DNAStringSet()
  
  
  # # # # # # # # # # # # # # # # # # # # 
  # read in PWMs (assumed jaspar)
  # # # # # # # # # # # # # # # # # # # # 
  print("reading PWMs")
  pfms_combined <- readJASPARMatrix(file_PWMs, matrixClass=c("PFM"))
  motif_names <- names(pfms_combined)
  PWMs_list <- c()
  for (motif in motif_names){
    pfm_oi <- pfms_combined[[motif]]
    PWMs_list <- c(PWMs_list,toPWM(pfm_oi))
  }
  PWMs <- do.call(PWMatrixList,PWMs_list)
  PWMs_ids <- ID(PWMs)
  names(PWMs) <- ID(PWMs)
  motif_color_map <- GLASBEYLUT$hex[1:(length(PWMs)+2)]
  motif_color_map <- motif_color_map[c(-1,-5)]
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # loop through sequences, find hits, print table and plot
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  if (!dir.exists("outs")){
    dir.create("outs")
  }
  
  
  print("working through sequences")
  
  df_PWM_scores <- data.frame()
  
  print(name_seqs)
  
  for (name_seq_oi in name_seqs){
    
    seq_oi <- seqs[name_seq_oi]
    
    df_PWM_scores_seq_oi <- data.frame()
    for (n_motif_oi in seq(length(PWMs))){
      PWM_oi <- PWMs[[n_motif_oi]]
      df_PWM_scores_oi <- PWM_scan_seq_oi_motif_oi_v2_20220206(seq_oi,PWM_oi)
      df_PWM_scores_oi$seq <- name_seq_oi
      df_PWM_scores_seq_oi <- rbind(df_PWM_scores_seq_oi,df_PWM_scores_oi) # will be slow if considering very large sequences
    }
    df_PWM_scores <- rbind(df_PWM_scores,df_PWM_scores_seq_oi)
    
    # options(dplyr.summarise.inform = FALSE)
    df_pos_hits <- get_position_motif_hit_df_no_offset_20220207(df_PWM_scores_seq_oi,name_seq_oi,score_threshold,PWMs,seq_oi)
    df_pos_hits_short <- df_pos_hits[[1]]
    df_pos_hits_long <- df_pos_hits[[2]]
    df_pos_hits_short$sequence <- name_seq_oi
    
     
     if (print_table_bool){
       date_str <- str_replace_all(Sys.Date(),"-","") 
       write.table(df_pos_hits_short,sprintf("outs/motif_hits_%s_score_thresh_%d_%s.txt",
                                             name_seq_oi,score_threshold,date_str), 
                   sep="\t", quote=FALSE, row.names=FALSE)
     }
    
    if (plot_bool){
      plt_seq_motif <- plot_sequence_with_motifs_20220202(seq_oi,df_pos_hits_long,motif_color_map)
      fig_name <- sprintf("outs/seq_w_motif_hits_%s_score_thresh_%d_%s.pdf",
                          name_seq_oi,score_threshold,date_str)
      pdf(fig_name,height=3,width=24)
      print(plt_seq_motif)
      dev.off()
    }
    

  }
  
  return(df_PWM_scores)
  
}