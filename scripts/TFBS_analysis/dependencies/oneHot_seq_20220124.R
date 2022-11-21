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