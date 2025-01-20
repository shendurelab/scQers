
# installation
library(ArchR)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(scales)
library(cowplot)
library(viridis)
library(tictoc)

path_dir <- "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC"
setwd(path_dir)



set.seed(1)
addArchRThreads(threads = 4)

# load genome
addArchRGenome("mm10")

inputFiles <- c("/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_1a_fragments.tsv.gz",
               "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_1b_fragments.tsv.gz",
               "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_2a_fragments.tsv.gz",
               "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_2b_fragments.tsv.gz")
names(inputFiles) <- c("fragments_1a","fragments_1b","fragments_2a","fragments_2b")



ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles <- c("fragments_1a.arrow","fragments_1b.arrow","fragments_2a.arrow","fragments_2b.arrow")

k_values <- c(3,10,30)
LSI_methods <- c(1,2,3)

# compute double score with default parameters
# doublet_scores_param_sweep <- matrix()
tic()
for (arrow_file in ArrowFiles){
  for (kval in k_values){
    for (LSI_m in LSI_methods){
      
      print(sprintf("doublet_enrichment_scATAC_%s_knn_LSI_kval_%g_LSImet_%g_20210210.txt",
                    substr(arrow_file,1,12),kval, LSI_m))
            
      doubScores <- addDoubletScores(
        input = arrow_file,
        k = kval, #Refers to how many cells near a "pseudo-doublet" to count.
        knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
        LSIMethod = LSI_m
      )
      
      print(sprintf("k-value: %g, LSI method: %g",
                    kval, LSI_m))
      print(toc())
      
      write.table(doubScores[[1]]$doubletEnrich,
                  sprintf("doublet_enrichment_scATAC_%s_knn_LSI_kval_%g_LSImet_%g_20210210.txt",
                          substr(arrow_file,1,12),kval, LSI_m), quote=FALSE, col.names=FALSE)
      
      # doublet_scores_param_sweep[kval, LSI_m, U_neigh, U_min] <- doubScores[[1]]$doubletEnrich
      
    }
  }
}




# no sweep for UMAP parameters this time around, found to not be critical on sweep from 02/10/2021
UMAP_neighbours <- 40 #c(10,40,100)
UMAP_min_dist <- 0.4 #c(0.1,0.4,1.2)

# compute double score with default parameters
# doublet_scores_param_sweep <- matrix()
tic()
for (arrow_file in ArrowFiles){
  for (kval in k_values){
    for (LSI_m in LSI_methods){
      # for (U_neigh in UMAP_neighbours){
      #   for (U_min in UMAP_min_dist){
      
      doubScores <- addDoubletScores(
        input = ArrowFiles[3],
        k = kval, #Refers to how many cells near a "pseudo-doublet" to count.
        knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
        LSIMethod = LSI_m,
        UMAPParams = list(n_neighbors = U_neigh, min_dist = U_min, metric = "euclidean", verbose =
                            FALSE)
      )
      
      print(sprintf("k-value: %g, LSI method: %g, UMAP neighbours: %g, UMAP min dist: %g",
                    kval, LSI_m, U_neigh, U_min))
      print(toc())
      
      write.table(doubScores[[1]]$doubletEnrich,
                  sprintf("doublet_enrichment_scATAC_%s_knn_UMAP_kval_%g_LSImet_%g_UMAPneighbours_%g_UMAPmindist_%g_20210210.txt",
                          arrow_file,kval, LSI_m, U_neigh, U_min), quote=FALSE, col.names=FALSE)
      

      #   }
      # }
    }
  }
}



# reading in data and comparing doublet scores

files <-  list.files(path = path_dir, pattern = glob2rx("doublet_enrichment_scATAC_EB1a*_20210204.txt"))
df1 <- data.frame(read.table(files[1],sep = " ", header = FALSE, comment.char = ""))

kval <- c()
LSImeth <- c()
UMAPneigh <- c()
UMAPmindist <- c()

df_doublet <- data.frame(Cell=df1$V1)
counter <- 1
for (file in files){

  df_oi <- data.frame(read.table(file,sep = " ", header = FALSE, comment.char = ""))
  colnames(df_oi) <- c("Cell",sprintf("DoubletEnrich%g",counter))
  
  df_doublet <- left_join(df_doublet,df_oi)
  
  
  # parsing params
  ind1 <- str_locate(file,"kval_")
  ind2 <- str_locate(file,"_LSImet_")
  ind3 <- str_locate(file,"_UMAPneighbours_")
  ind4 <- str_locate(file,"_UMAPmindist_")
  ind5 <- str_locate(file,"_20210204")
  
  kval[counter] <-     as.numeric(substring(file,ind1[1,2]+1,ind2[1,1]-1))
  LSImeth[counter] <-  as.numeric(substring(file,ind2[1,2]+1,ind3[1,1]-1))
  UMAPneigh[counter] <-  as.numeric(substring(file,ind3[1,2]+1,ind4[1,1]-1))
  UMAPmindist[counter] <-  as.numeric(substring(file,ind4[1,2]+1,ind5[1,1]-1))
  
  counter <- counter+1
  
}


files2 <-  list.files(path = path_dir, pattern = glob2rx("doublet_enrichment_scATAC_EB1a*_20210205.txt"))
df1 <- data.frame(read.table(files2[1],sep = " ", header = FALSE, comment.char = ""))

df_doublet2 <- data.frame(Cell=df1$V1)
counter <- 1
for (file in files2){
  
  df_oi <- data.frame(read.table(file,sep = " ", header = FALSE, comment.char = ""))
  colnames(df_oi) <- c("Cell",sprintf("DoubletEnrich_2_%g",counter))
  
  df_doublet2 <- left_join(df_doublet2,df_oi)
  
  
  # # parsing params
  # ind1 <- str_locate(file,"kval_")
  # ind2 <- str_locate(file,"_LSImet_")
  # ind3 <- str_locate(file,"_UMAPneighbours_")
  # ind4 <- str_locate(file,"_UMAPmindist_")
  # ind5 <- str_locate(file,"_20210204")
  # 
  # kval[counter] <-     as.numeric(substring(file,ind1[1,2]+1,ind2[1,1]-1))
  # LSImeth[counter] <-  as.numeric(substring(file,ind2[1,2]+1,ind3[1,1]-1))
  # UMAPneigh[counter] <-  as.numeric(substring(file,ind3[1,2]+1,ind4[1,1]-1))
  # UMAPmindist[counter] <-  as.numeric(substring(file,ind4[1,2]+1,ind5[1,1]-1))
  
  counter <- counter+1
  
}


df_doublet_combo <- left_join(df_doublet,df_doublet2)


# overlaying overlap doublet on distribution of doublet scores
xlims <- c(8e-3,5)


ind1 <- 8
ind2 <- 17

bool <- df_doublet_combo$Number.of.Overlaps>30
bool <- df_doublet_combo$p.value<0.01

length(which(bool))

plt_comp <- ggplot()+
  geom_hex(aes(x=df_doublet_combo[,ind1+1],y=df_doublet_combo[,ind2+1], color = ..count..), bins = 100)+
  geom_point(aes(x=df_doublet_combo[bool,ind1+1],y=df_doublet_combo[bool,ind2+1]),color = "red")+
  theme(legend.position="none")+
  scale_y_log10(breaks = xlabels,
                labels = trans_format("log10", math_format(10^.x)),
                limits = xlims) +
  scale_x_log10(breaks = xlabels,
                labels = trans_format("log10", math_format(10^.x)),
                limits = xlims) +
  coord_fixed() +
  labs(x = sprintf("doublet score 1"),
       y=  sprintf("doublet score 2"))
plt_comp



# see analysis_output_ATAC-DoubletDetector_20210205

cell_bc <- str_sub(df_doublet_combo$Cell,14,-1)

df_doublet_combo <- transform(df_doublet_combo,Barcode=cell_bc)

df_doublet_combo <- left_join(df_doublet_combo,df_ATAC_DD)


# pairwise correlation
R2 <- matrix(data = NA, 
             nrow = length(files)+length(files2),
             ncol = length(files)+ length(files2))

for (ind1 in seq(1,length(files)+ length(files2))){
  for (ind2 in seq(1,length(files)+ length(files2))){
    # if (ind1>ind2){
    x <- df_doublet_combo[,ind1+1]
    y <- df_doublet_combo[,ind2+1]
    bool <- x>0 & !is.na(x) & y>0 & !is.na(y)
    x2 <- log(x[bool])
    y2 <- log(y[bool])
    R2[ind1,ind2] <- cor(x2,y2)^2
    # }
  }
}

all_files <- c(files,files2)
rownames(R2) <- str_sub(all_files,32,-14)
colnames(R2) <- str_sub(all_files,32,-14)

require(gplots)

library(RColorBrewer)
colSide <- brewer.pal(9, "Set1")[R2_ArchR_DD]

# heatmap representation of matrix? 

library(scales)
trans <- div_gradient_pal(high="black",low="white", space="Lab")
cols <- trans(R2_ArchR_DD/max(R2_ArchR_DD))

trans2 <- div_gradient_pal(high="black",low="white", space="Lab")
cols_p <- trans2(R2_ArchR_DD_pval/max(R2_ArchR_DD_pval))

heatmap.2(R2,margins = c(10,15),trace="none", symm = TRUE, col=gray.colors(n=100,rev=TRUE),cexRow = 0.5, RowSideColors = cols, ColSideColors=cols_p)



# pairwise correlation
R2 <- matrix(data = NA, 
             nrow = length(files),
             ncol = length(files))

for (ind1 in seq(1,length(files))){
  for (ind2 in seq(1,length(files))){
    # if (ind1>ind2){
    x <- df_doublet[,ind1+1]
    y <- df_doublet[,ind2+1]
    bool <- x>0 & !is.na(x) & y>0 & !is.na(y)
    x2 <- log(x[bool])
    y2 <- log(y[bool])
    R2[ind1,ind2] <- cor(x2,y2)^2
    # }
  }
}

rownames(R2) <- str_sub(files,32,-14)
colnames(R2) <- str_sub(files,32,-14)

RowSideColors=colSide

require(gplots)
# heatmap representation of matrix? 
heatmap(R2,margins = c(10,10))
heatmap.2(R2,margins = c(30,35),trace="none", symm = TRUE, col=gray.colors(n=100,rev=TRUE),font=2)



R2_ArchR_DD <- matrix(data = NA, 
             nrow = length(all_files),
             ncol = 1)

R2_ArchR_DD_pval <- matrix(data = NA, 
                      nrow = length(all_files),
                      ncol = 1)

for (ind1 in seq(1,length(all_files))){
    # if (ind1>ind2){
    x <- df_doublet_combo[,ind1+1]
    
    y <- df_doublet_combo$Number.of.Overlaps+pseudo_c
    bool <- x>0 & !is.na(x) & y>0 & !is.na(y)
    x2 <- log(x[bool])
    y2 <- log(y[bool])
    R2_ArchR_DD[ind1] <- cor(x2,y2)^2
    
    
    y <- -log(df_doublet_combo$p.value)
    
    bool <- x>0 & !is.na(x) & y>0 & !is.na(y) & is.finite(y)
    x2 <- log(x[bool])
    y2 <- log(y[bool])
    R2_ArchR_DD_pval[ind1] <- cor(x2,y2)^2
    # }
  }





# comparing distributions

for (ind in seq(1,length(all_files))){

  xlims <- c(1e-3,10)
  # xlims <- c(0,3)
  
  xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1)
  n_bins2 <- 100
  plt_dist <- ggplot() +
    stat_bin(bins = n_bins2, aes(x=df_doublet_combo[,ind+1],y=..count..), geom="step")+
    scale_x_log10(breaks = xlabels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlims)+
    scale_y_log10()+
    annotation_logticks(sides = "lb")+
    # xlim(xlims)+
    labs(x="Doublet Enrichment",y="Cell counts")
  
  
  
  file_name <- str_replace(all_files[ind],".txt",".pdf")
  fig_name <- sprintf("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/analysis_doublet_param_sweep/log_%s",file_name)
  pdf(file = fig_name,width=10.5,height=7)
  print(plt_dist)
  dev.off()
  counter <- counter+1

}

for (ind in seq(1,length(all_files))){
  
  # xlims <- c(1e-3,10)
  xlims <- c(0,6)
  
  xplot <- seq(min(xlims),max(xlims),length.out = 500)
  CDF_doublet <- ecdf(df_doublet_combo[,ind+1])
  
  xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1)
  n_bins2 <- 100
  plt_CDF <- ggplot() +
    geom_step(aes(x=xplot,y=CDF_doublet(xplot)))+
    # scale_x_log10(breaks = xlabels, 
    #               labels = trans_format("log10", math_format(10^.x)), 
    #               limits = xlims)+
    # scale_y_log10()+
    # annotation_logticks(sides = "lb")+
    xlim(xlims)+
    labs(x="Doublet Enrichment",y="CDF (cell counts)")
  
  
  
  file_name <- str_replace(all_files[ind],".txt",".pdf")
  fig_name <- sprintf("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/analysis_doublet_param_sweep/linear_CDF_%s",file_name)
  pdf(file = fig_name,width=10.5,height=7)
  print(plt_CDF)
  dev.off()
  counter <- counter+1
  
}

# 
# CDF_frac <- ecdf(expression_phenotype$fraction_cells_exp)
# xlims_frac <- c(1e-4,1)
# xlabels_frac <- c(0.1,1,10,100,1000)*1E-4
# xplot_frac <- 10^(seq(-4,0,0.01))
# n_bins1 <- 50
# plt_frac <- ggplot() +
#   geom_step(aes(x=xplot_frac,y=CDF_frac(xplot_frac)))+
#   scale_x_log10(breaks = xlabels_frac, 
#                 labels = trans_format("log10", math_format(10^.x)), 
#                 limits = xlims_frac)+
#   annotation_logticks(sides = "b") + labs(x="Fraction of cells with detected expression",y="CDF Gene counts")





ind_ref <- 8
# k 10
# LSI method 1
# n neighbours 40
# min dist  0.4
bool <- df_doublet_combo$p.value<0.01

for (ind in seq(1,length(all_files))){
  
  xlims <- c(8e-3,5)
  # xlims <- c(0,3)
  xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1)

  plt_comp <- ggplot()+
      geom_hex(aes(x=df_doublet_combo[,ind_ref+1],y=df_doublet_combo[,ind+1], color = ..count..), bins = 100)+
      geom_point(aes(x=df_doublet_combo[bool,ind_ref+1],y=df_doublet_combo[bool,ind+1]),color = "red")+
    theme(legend.position="none")+
    scale_y_log10(breaks = xlabels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlims) +
    scale_x_log10(breaks = xlabels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlims) +
    coord_fixed() + annotation_logticks() +
    labs(x = sprintf("doublet score 1"),
         y=  sprintf("doublet score 2"))
  # 
  # 
  # xlims <- c(0,3)
  # # xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1)
  # 
  # plt_comp <- ggplot()+
  #   geom_hex(aes(x=df_doublet_combo[,ind_ref+1],y=df_doublet_comb[,ind+1], color = ..count..), bins = 100)+
  #   geom_point(aes(x=df_doublet_combo[bool,ind_ref+1],y=df_doublet_combo[bool,ind+1]),color = "red")+
  #   theme(legend.position="none")+
  #   xlim(xlims)+ylim(xlims)+
  #   # scale_y_log10(breaks = xlabels,
  #   #               labels = trans_format("log10", math_format(10^.x)),
  #   #               limits = xlims) +
  #   # scale_x_log10(breaks = xlabels,
  #   #               labels = trans_format("log10", math_format(10^.x)),
  #   #               limits = xlims) +
  #   coord_fixed() +
  #   labs(x = sprintf("doublet score 1"),
  #        y=  sprintf("doublet score 2"))
  
  
  # file_name <- str_replace(files[ind],"_20210204.txt","_20210205.pdf")
  fig_name <- sprintf("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/analysis_doublet_param_sweep/log_comparison_%g_vs_%g.pdf",ind,ind_ref)
  pdf(file = fig_name,width=10.5,height=7)
  print(plt_comp)
  dev.off()
  counter <- counter+1
  
}




pseudo_c <- 0.2

xlims <- c(8e-3,10)
# xlims <- c(0,3)

xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)

ylims <- c(0.1,5e2)


for (ind in seq(1,length(all_files))){
  
plt_comp <- ggplot()+
  geom_hex(aes(x=df_doublet_combo[,ind+1],y=df_doublet_combo$Number.of.Overlaps+pseudo_c, color = ..count..), bins = 100)+
  theme(legend.position="none")+
  scale_y_log10(breaks = xlabels,
                labels = trans_format("log10", math_format(10^.x)),
                limits = ylims) +
  scale_x_log10(breaks = xlabels,
                labels = trans_format("log10", math_format(10^.x)),
                limits = xlims) +
  coord_fixed() + annotation_logticks() +
  labs(x = sprintf("doublet score ArchR"),
       y=  sprintf("doublet score: ATAC-DoubletDetector"))+
  labs(title=all_files[ind])

  file_name <- str_replace(all_files[ind],".txt",".pdf")
  
  plt_comp
  bla

  fig_name <- sprintf("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/analysis_doublet_param_sweep/comparison_ArchR_DD_%s.pdf",file_name)
  pdf(file = fig_name,width=10.5,height=7)
  print(plt_comp)
  dev.off()
  
  counter <- counter+1


}




pseudo_c <- 0.2

xlims <- c(8e-3,5)
# xlims <- c(0,3)

xlabels <- c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)
ylims <- c(0.1,5e2)

for (ind in seq(1,length(all_files))){
  
  plt_comp <- ggplot()+
    geom_hex(aes(x=df_doublet_combo[,ind+1],y=-log(df_doublet_combo$p.value), color = ..count..), bins = 100)+
    theme(legend.position="none")+
    scale_y_log10(breaks = xlabels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = ylims) +
    scale_x_log10(breaks = xlabels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlims) +
    coord_fixed() + annotation_logticks() +
    labs(x = sprintf("doublet score ArchR"),
         y=  sprintf("doublet score: ATAC-DoubletDetector -log p-value"))+
    labs(title=all_files[ind])
  
  file_name <- str_replace(all_files[ind],".txt",".pdf")
  
  fig_name <- sprintf("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/analysis_doublet_param_sweep/comparison_ArchR_DDpval_%s.pdf",file_name)
  pdf(file = fig_name,width=10.5,height=7)
  print(plt_comp)
  dev.off()
  
  counter <- counter+1
  
}










  # creating the main ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "SR_EB_all",
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)

# from the project, extract the doublet enrichment score.
# plot distribution for each sample. 


doublet_enrich <- proj$DoubletEnrichment
cell_samples <- proj$Sample

samples <- c("fragments_1a","fragments_1b","fragments_2a","fragments_2b")

n_cell_cut <- c(2418,2367,1317,1512)
names(n_cell_cut) <- samples

for (sample_oi in samples){
  # sample_oi <- "fragments_2b"
  id_sample_oi <- which(cell_samples==sample_oi)
  sum(doublet_enrich[id_sample_oi]>1)
  
  n_cells <- length(id_sample_oi)
  n_cells
  
  # cutoff assuming from Poisson expected doublets
  n_beads_tot <- 1E5
  lambda <- -log(1-n_cells/n_beads_tot)
  n_expected_doublets <- (1-exp(-lambda) - lambda*exp(-lambda))/(lambda*exp(-lambda))*n_cells
  n_expected_doublets
  
  # identifying enrichment cutoff:
  thresh <- seq(from = 0,to = 3,by = 0.001)
  n_cells_above <- vector(mode = "numeric", length = length(thresh))
  for (ind in seq(1,length(thresh))){
    n_cells_above[ind] = sum(doublet_enrich[id_sample_oi]>thresh[ind])
  }
  thresh_doublet <- thresh[max(which(n_cells_above>n_expected_doublets))]
  thresh_doublet_ArchR <- thresh[max(which(n_cells_above>n_cell_cut[sample_oi]))]
  
  
  thresh_doublet
  
  xlims <- c(0,3)
  ylims_labs <- c(1, 10, 100, 1000)
  ylims_hist <- c(1,2000)
  plt2 <- ggplot() + 
    geom_line(aes(x = c(thresh_doublet,thresh_doublet), y = ylims_hist), colour = "red", linetype = "dashed")+
    geom_line(aes(x = c(thresh_doublet_ArchR,thresh_doublet_ArchR), y = ylims_hist), colour = "blue", linetype = "dotted")+
    stat_bin(bins = 100, aes(x=doublet_enrich[id_sample_oi],y=..count..), geom="step")+
    scale_y_log10(breaks = ylims_labs,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = ylims_hist) +
    annotation_logticks(sides="l")+
    xlim(xlims)+
    labs(x = "Doublet enrichment score", y = "Cell counts", title = sprintf("%s, n = %d, log", sample_oi, n_cells))
  
  plt1 <- ggplot() + 
    geom_line(aes(x = c(thresh_doublet,thresh_doublet), y = ylims_hist), colour = "red", linetype = "dashed")+
    geom_line(aes(x = c(thresh_doublet_ArchR,thresh_doublet_ArchR), y = ylims_hist), colour = "blue", linetype = "dotted")+
    stat_bin(bins = 100, aes(x=doublet_enrich[id_sample_oi],y=..count..), geom="step")+
    xlim(xlims)+
    # annotation_logticks(sides="bl")+
    labs(x = "Doublet enrichment score", y = "Cell counts", title = sprintf("%s, n = %d",sample_oi, n_cells))
  
  plt <- plot_grid(plt1,plt2)
  
  
  fig_name <- sprintf("double_threshold_%s_20201130.pdf",sample_oi)
  pdf(file = fig_name,width=8,height=3)
  print(plt)
  dev.off()
}

# basic interaction with the project to pull information from memory
proj_nodoub <- filterDoublets(proj)
proj_nodoub <- saveArchRProject(ArchRProj = proj_nodoub, outputDirectory = "Save-ProjNoDoub", load = TRUE)


proj_nodoub <- loadArchRProject("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/Save-ProjNoDoub")

idxPass <- which(proj_nodoub$TSSEnrichment >= 8 & log10(proj_nodoub$nFrags)>3.3)

cellsPass <- proj_nodoub$cellNames[idxPass]

proj_nodoub_filter <- proj_nodoub[cellsPass, ]



proj_nodoub_filter <- addIterativeLSI(
  ArchRProj = proj_nodoub_filter,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:100,
  force = TRUE
)



LSI_mat <- getReducedDims(
  ArchRProj = proj_nodoub_filter,
  reducedDims = "IterativeLSI",
  returnMatrix = TRUE,
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75
)


proj_nodoub_filter <- addClusters(
  input = proj_nodoub_filter,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)



cM <- confusionMatrix(paste0(proj_nodoub_filter$Clusters), paste0(proj_nodoub_filter$Sample))






geneScoreMat <- getMatrixFromArrow(
  ArrowFile = ArrowFiles,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  cellNames = NULL,
  ArchRProj = NULL,
  verbose = TRUE,
  binarize = FALSE,
  logFile = createLogFile("getMatrixFromArrow")
)





