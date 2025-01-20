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
require(gplots)
library(scales)
library(RColorBrewer)

path_dir <- "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC"
setwd(path_dir)

devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.1", repos = BiocManager::repositories())


set.seed(1)
addArchRThreads(threads = 4)

# load genome
addArchRGenome("mm10")

# # generated around 02/16/2021
# inputFiles <- c("/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_1a_fragments.tsv.gz",
#                 "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_1b_fragments.tsv.gz",
#                 "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_2a_fragments.tsv.gz",
#                 "/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/sample_2b_fragments.tsv.gz")
# names(inputFiles) <- c("fragments_1a","fragments_1b","fragments_2a","fragments_2b")
# 
# 
# 
# ArrowFiles <- createArrowFiles(
#   inputFiles = inputFiles,
#   sampleNames = names(inputFiles),
#   minTSS = 4, #Dont set this too high because you can always increase later
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )


ArrowFiles <- c("fragments_1a.arrow","fragments_1b.arrow","fragments_2a.arrow","fragments_2b.arrow")

# creating the main ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "SR_EB_all_v2_analysis",
  copyArrows = FALSE # TRUE is recommended so that you maintain an unaltered copy for later usage.
)


# this is NOT what I want, because it is not the final parameters that I want to use, it was from the parameter sweep that I did 
doublet_enrich <- proj$DoubletEnrichment

# # these are the best parameters, as determined by correlation to the DoubletDetector 
# kval <- 10
# LSI_m <- 1
# knnMethod <- "LSI"
# doubScores <- addDoubletScores(
#   input = arrow_file,
#   k = kval, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = knnMethod, #Refers to the embedding to use for nearest neighbor search.
#   LSIMethod = LSI_m # normalization of the LSI matrix 
# )


cell_samples <- proj$Sample


# basic QC
idxPass <- which(proj$TSSEnrichment >= 8 & log10(proj$nFrags)>3.3)
cellsPass <- proj$cellNames[idxPass]
proj_TSS_nfrag_QC <- proj[cellsPass, ]


# dimensionality reduction: robustness there? 
proj_TSS_nfrag_QC <- addIterativeLSI(
  ArchRProj = proj_TSS_nfrag_QC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:30,
  force = TRUE  # overwritting ArchR project
)




# some incompatbility with the latest uwot version. Trying to install older one...
# c.f., https://github.com/GreenleafLab/ArchR/issues/424
# think it is resolved now
# devtools::install_version("RcppAnnoy", version = "0.0.16") #needed to remove/rename my file ~/.R/Makevars for this to work
# devtools::install_version("uwot", version = "0.1.8") #Note - when prompted to update packages, answer 3 (None) - thus avoiding replacing RcppAnnoy v0.0.16 with v0.0.17


# starting point for 03/16/2021
# proj_TSS_nfrag_QC <- saveArchRProject(ArchRProj = proj_TSS_nfrag_QC, outputDirectory = "SR_EB_all_v2_analysis", load = TRUE)
proj_TSS_nfrag_QC <- loadArchRProject("/Users/jbl/Documents/UW/projects/EB_sc_data/scATAC/SR_EB_all_v2_analysis")

proj_TSS_nfrag_QC <- addClusters(
    input = proj_TSS_nfrag_QC,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.5,
    force=TRUE)

# low dimensional UMAP
proj_TSS_nfrag_QC <- addUMAP(
  ArchRProj = proj_TSS_nfrag_QC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


# plotting
p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

head(proj_TSS_nfrag_QC$cellNames)

# add the doublet score information to the UMAP. How? 
files_doublet_enrichment_scores_ArchR <- c("doublet_enrichment_scATAC_fragments_1a_knn_LSI_kval_10_LSImet_1_20210210.txt",
                                     "doublet_enrichment_scATAC_fragments_1b_knn_LSI_kval_10_LSImet_1_20210210.txt",
                                     "doublet_enrichment_scATAC_fragments_2a_knn_LSI_kval_10_LSImet_1_20210210.txt",
                                     "doublet_enrichment_scATAC_fragments_2b_knn_LSI_kval_10_LSImet_1_20210210.txt")
# reading the data
df_doublet_score <- c()
for (file in files_doublet_enrichment_scores_ArchR){
  df_doublet_score_oi <-  read.table(file,sep=" ",comment.char = "")
  
  df_doublet_score <- rbind(df_doublet_score,df_doublet_score_oi)
}

bool_cells_retained <- df_doublet_score$V1 %in% proj_TSS_nfrag_QC$cellNames

proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = log10(df_doublet_score$V2[bool_cells_retained]),
  name = "ArchR_doublet_knn_LSI_kval_10_LSImet_1",
  cells = df_doublet_score$V1[bool_cells_retained],
  force = TRUE
)

# add the doublet score information to the UMAP. How? 
files_doublet_enrichment_scores_ArchR <- c("doublet_enrichment_scATAC_fragments_1a_knn_UMAP_kval_10_LSImet_1_UMAPneighbours_40_UMAPmindist_0.4_20210216.txt",
                                           "doublet_enrichment_scATAC_fragments_1b_knn_UMAP_kval_10_LSImet_1_UMAPneighbours_40_UMAPmindist_0.4_20210216.txt",
                                           "doublet_enrichment_scATAC_fragments_2a_knn_UMAP_kval_10_LSImet_1_UMAPneighbours_40_UMAPmindist_0.4_20210216.txt",
                                           "doublet_enrichment_scATAC_fragments_2b_knn_UMAP_kval_10_LSImet_1_UMAPneighbours_40_UMAPmindist_0.4_20210216.txt")
# reading the data
df_doublet_score <- c()
for (file in files_doublet_enrichment_scores_ArchR){
  df_doublet_score_oi <-  read.table(file,sep=" ",comment.char = "")
  
  df_doublet_score <- rbind(df_doublet_score,df_doublet_score_oi)
}
bool_cells_retained <- df_doublet_score$V1 %in% proj_TSS_nfrag_QC$cellNames
proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = log10(df_doublet_score$V2[bool_cells_retained]),
  name = "ArchR_doublet_knn_UMAP_kval_10_LSImet_1",
  cells = df_doublet_score$V1[bool_cells_retained],
  force = TRUE
)

plt_comp <- ggplot()+
  geom_hex(aes(x=proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1,y=proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1, color = ..count..), bins = 100)+
  theme(legend.position="none")+
  coord_fixed()
plt_comp


# plotting
p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", 
                    name = "ArchR_doublet_knn_LSI_kval_10_LSImet_1", embedding = "UMAP",plotAs="points",colorLimits=c(-2,0))
ggAlignPlots(p1, p2, type = "h")


p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", 
                    name = "ArchR_doublet_knn_UMAP_kval_10_LSImet_1", embedding = "UMAP",plotAs="points")
ggAlignPlots(p1, p2, type = "h")


sum(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1>0)/length(proj_TSS_nfrag_QC$cellNames)
sum(proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1>0)/length(proj_TSS_nfrag_QC$cellNames)


bool_LSI_high_UMAP_low <- as.numeric(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1>0 & proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1<0)
bool_LSI_low_UMAP_high <- as.numeric(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1<0 & proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1>0)
bool_LSI_high_UMAP_high <- as.numeric(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1>0 & proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1>0)
bool_LSI_low_UMAP_low <- as.numeric(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1<0 & proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1<0)

proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = bool_LSI_high_UMAP_low,
  name = "bool_LSI_high_UMAP_low",
  cells = proj_TSS_nfrag_QC$cellNames,
  force = TRUE
)
proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = bool_LSI_low_UMAP_high,
  name = "bool_LSI_low_UMAP_high",
  cells = proj_TSS_nfrag_QC$cellNames,
  force = TRUE
)
proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = bool_LSI_high_UMAP_high,
  name = "bool_LSI_high_UMAP_high",
  cells = proj_TSS_nfrag_QC$cellNames,
  force = TRUE
)
proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = bool_LSI_low_UMAP_low,
  name = "bool_LSI_low_UMAP_low",
  cells = proj_TSS_nfrag_QC$cellNames,
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "bool_LSI_high_UMAP_low", embedding = "UMAP",randomize=FALSE,plotAs="points")
p2 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "bool_LSI_low_UMAP_high", embedding = "UMAP",randomize=FALSE,plotAs="points")
p3 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "bool_LSI_high_UMAP_high", embedding = "UMAP",randomize=FALSE,plotAs="points")
p4 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "bool_LSI_low_UMAP_low", embedding = "UMAP",randomize=FALSE,plotAs="points")

ggAlignPlots(p1, p2,p3, p4, type = "h")

proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = log10(proj_TSS_nfrag_QC$nFrags),
  name = "log10_nFrags",
  cells=proj_TSS_nfrag_QC$cellNames,
  force = TRUE
)

plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "log10_nFrags", embedding = "UMAP",randomize=FALSE, plotAs="points",colorLimits=c(3.3,4))


# adding doublet detector values
files_DoubletDector_scores <- c("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/DoubletDetector_outputs/out_scATAC_1a_20210121/DoubletProbabilities.txt",
                                "/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/DoubletDetector_outputs/out_scATAC_1b_20210207/DoubletProbabilities.txt",
                                "/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/DoubletDetector_outputs/out_scATAC_2a_20210207/DoubletProbabilities.txt",
                                "/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/DoubletDetector_outputs/out_scATAC_2b_20210207/DoubletProbabilities.txt")
sample_names <- c("fragments_1a","fragments_1b","fragments_2a","fragments_2b")
names(sample_names) <-files_DoubletDector_scores

id_file <- 1
df_doublet_detector_oi <-  read.table(files_DoubletDector_scores[id_file],sep="\t",comment.char = "",header=TRUE)
cell_names <- paste0(sample_names[files_DoubletDector_scores[id_file]],"#",df_doublet_detector_oi$barcode)
sum(cell_names %in% proj_TSS_nfrag_QC$cellNames)/length(cell_names)

df_doublet_detector <- c()
for (file in files_DoubletDector_scores){
  df_doublet_detector_oi <-  read.table(file,sep="\t",comment.char = "",header=TRUE)
  cell_names <- paste0(sample_names[file],"#",df_doublet_detector_oi$barcode)
  
  bool_cells <- cell_names %in% proj_TSS_nfrag_QC$cellNames
  print(sum(bool_cells))

  df_doublet_detector_oi <- transform(df_doublet_detector_oi,Cell_names=cell_names)
  
  df_doublet_detector <- rbind(df_doublet_detector,df_doublet_detector_oi)
}

# adding the doublet scores as cellColData, just need to name cells properly from the doublet data.
# what are cell names in current project file?

# basic interaction with the project to pull information from memory
# proj_nodoub <- filterDoublets(proj)
# proj_nodoub <- saveArchRProject(ArchRProj = proj_nodoub, outputDirectory = "Save-ProjNoDoub", load = TRUE)
# proj_nodoub <- loadArchRProject("/Users/jbl/Documents/UW/projects/scATAC_analysis/SR_EB_10x_scATAC/Save-ProjNoDoub")

bool_cells_retained <- df_doublet_detector$Cell_names %in% proj_TSS_nfrag_QC$cellNames

sum(df_doublet_detector$Cell_names %in% proj_TSS_nfrag_QC$cellNames)/length(df_doublet_detector$Cell_names)
sum(df_doublet_detector$Cell_names %in% proj_TSS_nfrag_QC$cellNames)/length(proj_TSS_nfrag_QC$cellNames)

proj_TSS_nfrag_QC <- addCellColData(
  ArchRProj = proj_TSS_nfrag_QC,
  data = log10(-log10(df_doublet_detector$p.value[bool_cells_retained])),
  name = "DoubletDector_log10_minus_log10_p_val",
  cells = df_doublet_detector$Cell_names[bool_cells_retained],
  force = TRUE
)

min(proj_TSS_nfrag_QC$DoubletDector_log10_minus_log10_p_val,na.rm=TRUE)

proj_TSS_nfrag_QC$DoubletDector_log10_minus_log10_p_val[is.na(proj_TSS_nfrag_QC$DoubletDector_log10_minus_log10_p_val)] <- -1

p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", 
                    name = "DoubletDector_log10_minus_log10_p_val", embedding = "UMAP",plotAs="points",colorLimits=c(-0.7,0.5))
ggAlignPlots(p1, p2, type = "h")


sum(is.na(proj_TSS_nfrag_QC$DoubletDector_minus_log10_p_val))
bool_DD <- is.na(proj_TSS_nfrag_QC$DoubletDector_minus_log10_p_val)

quantile(proj_TSS_nfrag_QC$nFrags[bool_DD])
quantile(proj_TSS_nfrag_QC$nFrags[!bool_DD])

plt_comp <- ggplot()+
  geom_hex(aes(x=proj_TSS_nfrag_QC$log10_nFrags,y=proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1, color = ..count..), bins = 100)+
  theme(legend.position="none")+
  coord_fixed()
plt_comp

plt_comp <- ggplot()+
  geom_hex(aes(x=proj_TSS_nfrag_QC$log10_nFrags,y=proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1, color = ..count..), bins = 100)+
  theme(legend.position="none")
plt_comp

plt_comp <- ggplot()+
  geom_hex(aes(x=proj_TSS_nfrag_QC$log10_nFrags,y=log10(proj_TSS_nfrag_QC$DoubletDector_minus_log10_p_val), color = ..count..), bins = 100)+
  theme(legend.position="none")
plt_comp



# more granular view with stratificatio by cluster. Violin plot etc? 
p2 <- plotGroups(
  ArchRProj = proj_TSS_nfrag_QC, 
  groupBy = "Clusters", 
  colorBy = "cellColData", 
  name = "ArchR_doublet_knn_UMAP_kval_10_LSImet_1",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj_TSS_nfrag_QC, 
  groupBy = "Clusters", 
  colorBy = "cellColData", 
  name = "ArchR_doublet_knn_LSI_kval_10_LSImet_1",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4 <- plotGroups(
  ArchRProj = proj_TSS_nfrag_QC, 
  groupBy = "Clusters", 
  colorBy = "cellColData", 
  name = "DoubletDector_log10_minus_log10_p_val",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p1 <- plotGroups(
  ArchRProj = proj_TSS_nfrag_QC, 
  groupBy = "Clusters", 
  colorBy = "cellColData", 
  name = "log10_nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
ggAlignPlots(p1,p2, p3,p4, type = "v")


p2
p1 <- plotEmbedding(ArchRProj = proj_TSS_nfrag_QC, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")


# putative doublet clusters (with resolution=0.5)
bad_clusters <- c("C3","C8","C5","C6","C7")
bool_cluster <- proj_TSS_nfrag_QC$Clusters %in% bad_clusters
sum(bool_cluster)/length(bool_cluster)

thresh_log10_nFrags <- 4.25
bool_nFrags <- proj_TSS_nfrag_QC$log10_nFrags>thresh_log10_nFrags
sum(bool_nFrags)/length(bool_nFrags)

LSI_doublet_thresh <- 0
bool_LSI_doublet <- proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1>LSI_doublet_thresh
sum(bool_LSI_doublet)/length(bool_LSI_doublet)

UMAP_doublet_thresh <- 0
bool_UMAP_doublet <- proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1>UMAP_doublet_thresh
sum(bool_UMAP_doublet)/length(bool_UMAP_doublet)


DoubletDetector_thresh <- 0.3
bool_DoubletDetector <- proj_TSS_nfrag_QC$DoubletDector_log10_minus_log10_p_val>DoubletDetector_thresh
sum(bool_DoubletDetector)/length(bool_DoubletDetector)

bool_individual_thresh <- bool_DoubletDetector | bool_UMAP_doublet | bool_LSI_doublet | bool_nFrags
sum(bool_individual_thresh)/length(bool_individual_thresh)


bool_cluster_individual_thresh <- bool_individual_thresh | bool_cluster
sum(bool_cluster_individual_thresh)/length(bool_cluster_individual_thresh)


quantile(proj_TSS_nfrag_QC$nFrags,seq(0,1,0.05))
quantile(proj_TSS_nfrag_QC$ArchR_doublet_knn_UMAP_kval_10_LSImet_1,seq(0,1,0.05))
quantile(proj_TSS_nfrag_QC$ArchR_doublet_knn_LSI_kval_10_LSImet_1,seq(0,1,0.05))
quantile(proj_TSS_nfrag_QC$DoubletDector_log10_minus_log10_p_val,seq(0,1,0.05))



## removing putative doublets
bool_cluster_individual_thresh

# saved at the end here. 
proj_TSS_nfrag_QC_no_doub <- proj_TSS_nfrag_QC[which(!bool_cluster_individual_thresh), ]
proj_TSS_nfrag_QC_no_doub <- saveArchRProject(ArchRProj = proj_TSS_nfrag_QC_no_doub, outputDirectory = "SR_EB_all_samples_doublet_removed", load = TRUE)



# clusters in the UMAP vs. sample rep identity? tech. vs. biol. rep?
cM <- confusionMatrix(paste0(proj_TSS_nfrag_QC$Clusters), paste0(proj_TSS_nfrag_QC$Sample))


library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)






