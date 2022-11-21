PWM_scan_seq_oi_motif_oi_v2_20220206 <- function(seq_oi,PWM_oi){
  
  # get reverse complement and convert to one-hot encoding
  seq_oi_RC <- reverseComplement(seq_oi)
  seq_size <- width(seq_oi)
  oh_seq <- Reshape(oneHot_seq_20220124(seq_oi),4,seq_size)
  oh_seq_RC <- Reshape(oneHot_seq_20220124(seq_oi_RC),4,seq_size)
  
  # PWM_oi <- PWMs_combined[["GATA4_MOUSE.H11MO.0.A"]]
  # PWM_oi <- PWMs_combined[["Sox17_MA0078.1"]]
  
  # PWM_oi <- Matrix(motif_oi)/colSums(motif_oi)
  name_motif <- name(PWM_oi)
  PWM_oi <- Matrix(PWM_oi)
  
  region_size <- length(oh_seq[1,])

  df_PWM_scores <- data.frame()
  
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
  
  df_PWM_oi_for <- data.frame(relative_positions=1:length(convo), score=convo, strand="+", motif=name_motif)
  df_PWM_oi_rev <- data.frame(relative_positions=1:length(convo), score=convo_reverse, strand="-", motif=name_motif)
  df_PWM_oi_max <- data.frame(relative_positions=1:length(convo), score=pmax(convo,convo_reverse), strand="max", motif=name_motif)
  
  df_PWM_scores <- rbind(df_PWM_oi_for,df_PWM_oi_rev,df_PWM_oi_max)
  
  return(df_PWM_scores)

  
}
  