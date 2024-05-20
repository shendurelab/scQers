library(Seurat)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(scales)
library(tictoc)
library(irlba)
# library(schex) # might need to install

library(DropletUtils)

library(pryr)


setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/")

wd_oi <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608"

# read raw data
reps <- c("2B1","2B2")



# # # PRE FILTERING
mt_thresh <- c(1,15)
UMI_thresh <- c(1000,400)
names(UMI_thresh) <- reps

for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_%s",rep_oi)
  gex_scREA_mat <- Read10X(data.dir =path_10x_mtx)
  
  
  # create separate Seurat objects
  gex_scREA <- CreateSeuratObject(counts = gex_scREA_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex_scREA[["percent_mt"]] <- PercentageFeatureSet(gex_scREA, pattern = "^mt-")
  
  ## QC plot
  # fig_QC <- ggplot()+geom_hex(aes(x=gex_scREA$nCount_RNA,y=gex_scREA$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+scale_x_log10()+scale_y_log10()+
  #   geom_vline(xintercept=UMI_thresh[rep_oi], color="red", linetype="dashed")+
  #   geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")
  # 
  #   
  # fig_QC2 <- ggplot()+geom_step(aes(x=gex_scREA$nCount_RNA,y=1-..y..),stat="ecdf")+
  #   scale_x_log10()+scale_y_log10()+
  #   geom_vline(xintercept=UMI_thresh[rep_oi], color="red", linetype="dashed")
  # 
  # 
  # 
  # fig_QC_name <- sprintf("mt_vs_RNA_UMI_gex_scREA_%s_20220610.pdf",rep_oi)
  # 
  # pdf(fig_QC_name,width=8,height=4)
  # print(fig_QC+fig_QC2)
  # dev.off()
  
  
  # subsetting on cells
  gex_scREA_2 <- subset(gex_scREA, subset = (percent_mt < mt_thresh[2]) & (percent_mt > mt_thresh[1]) & (nCount_RNA >  UMI_thresh[rep_oi] ) )
  
  print(dim(gex_scREA))
  print(dim(gex_scREA_2))
  
  write10xCounts(sprintf("%s/gex_%s_filtered_10x_mtx_20220610/",wd_oi,rep_oi),
                 gex_scREA_2@assays$RNA@counts,
                 barcodes = colnames(gex_scREA_2),
                 gene.id = rownames(gex_scREA_2))
  
}

#gex_scREA_2 <- subset(gex_scREA, subset = (percent_mt < mt_thresh[2]) & (nCount_RNA >  UMI_thresh ) )



# run scrublet, see scrublet_scREA_seq031_20220610.py

# read raw data
reps <- c("A1","A2","B1","B2","2B1","2B2")


files_data <- c("/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_A1_filtered_10x_mtx_20220525",
                "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_A2_filtered_10x_mtx_20220525",
                "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_B1_filtered_10x_mtx_20220525",
                "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_B2_filtered_10x_mtx_20220525",
                "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_2B1_filtered_10x_mtx_20220610",
                "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_2B2_filtered_10x_mtx_20220610")
names(files_data) <- reps

files_doublet <- c("/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_A1_filtered_10x_mtx_20220525/scrublet_score_scRNAgex_A1_filtered_10x_mtx_20220525.txt",
                   "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_A2_filtered_10x_mtx_20220525/scrublet_score_scRNAgex_A2_filtered_10x_mtx_20220525.txt",
                   "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_B1_filtered_10x_mtx_20220525/scrublet_score_scRNAgex_B1_filtered_10x_mtx_20220525.txt",
                   "/Users/jbl/Documents/UW/Data/sequencing_runs/seq028_mEB_sc_rep_v1_20220524/gex_outs_20220525/gex_B2_filtered_10x_mtx_20220525/scrublet_score_scRNAgex_B2_filtered_10x_mtx_20220525.txt",
                   "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_2B1_filtered_10x_mtx_20220610/scrublet_score_scRNAgex_2B1_filtered_10x_mtx_20220610.txt",
                   "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_2B2_filtered_10x_mtx_20220610/scrublet_score_scRNAgex_2B2_filtered_10x_mtx_20220610.txt")
names(files_doublet) <- reps

# # # READ ALL W/ DOUBLET SCORES
for (rep_oi in reps){
  
  print(rep_oi)
  str_file <- files_data[rep_oi] %>% as.character
  gex_scREA_mat <- Read10X(data.dir = str_file)
  
  # create separate Seurat objects
  gex_scREA <- CreateSeuratObject(counts = gex_scREA_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex_scREA[["percent_mt"]] <- PercentageFeatureSet(gex_scREA, pattern = "^mt-")
  
  df_doublet <- read.table(files_doublet[rep_oi],sep="\t",header=TRUE)
  
  rownames(df_doublet) <- df_doublet$cell_bc
  gex_scREA$doublet_score <- df_doublet[Cells(gex_scREA),2]
  
  assign(sprintf("gex_scREA_%s",rep_oi),gex_scREA)

}

# # # # threshold on scrublet score
scrublet_thresh <- 0.3
gex_scREA_A1_2 <- subset(gex_scREA_A1, subset = doublet_score<scrublet_thresh) 
gex_scREA_A2_2 <- subset(gex_scREA_A2, subset = doublet_score<scrublet_thresh) 
gex_scREA_B1_2 <- subset(gex_scREA_B1, subset = doublet_score<scrublet_thresh) 
gex_scREA_B2_2 <- subset(gex_scREA_B2, subset = doublet_score<scrublet_thresh) 

gex_scREA_2B1_2 <- subset(gex_scREA_2B1, subset = doublet_score<scrublet_thresh) 
gex_scREA_2B2_2 <- subset(gex_scREA_2B2, subset = doublet_score<scrublet_thresh) 


# # # MERGE ALL IN ONE SEURAT OBJECT

Idents(gex_scREA_A1_2) <- "A1"
Idents(gex_scREA_A2_2) <- "A2"
Idents(gex_scREA_B1_2) <- "B1"
Idents(gex_scREA_B2_2) <- "B2"

Idents(gex_scREA_2B1_2) <- "2B1"
Idents(gex_scREA_2B2_2) <- "2B2"


gex_scREA_all <- merge(gex_scREA_A1_2, y = c(gex_scREA_A2_2,
                                             gex_scREA_B1_2,
                                             gex_scREA_B2_2,
                                             gex_scREA_2B1_2,
                                             gex_scREA_2B2_2),
                       add.cell.ids = c("A1","A2","B1","B2","2B1", "2B2"))



# normalization/dimensional reduction
gex_scREA_3 <- gex_scREA_all
gex_scREA_3 <- NormalizeData(gex_scREA_3, normalization.method = "LogNormalize", scale.factor = 10000)
gex_scREA_3 <- FindVariableFeatures(gex_scREA_3, selection.method = "vst", nfeatures = 1000, verbose = TRUE)
all.genes <- rownames(gex_scREA_3)
gex_scREA_3 <- ScaleData(gex_scREA_3, features = all.genes)
gex_scREA_3 <- RunPCA(gex_scREA_3, features = VariableFeatures(object = gex_scREA_3), verbose = FALSE, npcs = 100)



gex_scREA_4 <- gex_scREA_3
top_pc <- 50 # before was 40
gex_scREA_4 <- FindNeighbors(gex_scREA_4, dims = 1:top_pc)
gex_scREA_4 <- FindClusters(gex_scREA_4, resolution = 0.2) # before was 0.5
gex_scREA_4 <- RunUMAP(gex_scREA_4, dims = 1:top_pc, n.neighbors = 50, seed.use = 42)

plt0 <- DimPlot(gex_scREA_4, reduction = "umap", label=TRUE)
plt1 <- FeaturePlot(object = gex_scREA_4,  features = 'nCount_RNA',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
plt2 <- DimPlot(object = gex_scREA_4,  split.by='orig.ident',group.by='orig.ident')
plt0+plt1+plt2+plot_layout(design="AABB##
                                   AABB##
                                   CCCCCC")


DimPlot(object = gex_scREA_4,  split.by = 'orig.ident')
FeaturePlot(object = gex_scREA_4,  slot = "count", features = c("Lama1","Sox17","Gata4","Foxa2","Gsc","T","Eomes","Tubb2b","Col5a1","Col1a1","Sox2","Nanog","Pou5f1"),  pt.size = 0.02,  max.cutoff = 'q95',  order = TRUE)
FeaturePlot(object = gex_scREA_4,  slot = "count", features = c("Lama1","Sox17","Gata4","Foxa2","Tubb2b","Col5a1","Col1a1","Sox2","Nanog","Pou5f1"),  pt.size = 0.02,  max.cutoff = 'q95',  order = TRUE)


# cluster representation:
table(gex_scREA_4$seurat_clusters,gex_scREA_4$orig.ident)/colSums(table(gex_scREA_4$seurat_clusters,gex_scREA_4$orig.ident))
table(gex_scREA_4$seurat_clusters,gex_scREA_4$orig.ident)



# saving some outs
saveRDS(gex_scREA_4,"gex_scREA_mEB_bottlenecked_all_reps_A_B_2B_20220610.RDS")

gex_scREA_4 <- readRDS("/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/gex_scREA_mEB_bottlenecked_all_reps_A_B_2B_20220610.RDS")

raw_cellBC_pdT <- str_split(Cells(gex_scREA_4),"_") %>% lapply("[[",2) %>% substr(1,16)
df_metadata_cells <- data.frame(cellBC=Cells(gex_scREA_4),
                                raw_cellBC_pdT,
                                rep_id=gex_scREA_4$orig.ident,
                                gex_UMI=gex_scREA_4$nCount_RNA,
                                gex_percent_mt=round(gex_scREA_4$percent_mt,2),
                                cluster_id=gex_scREA_4$seurat_clusters)

write.table(df_metadata_cells,"metadata_filtered_cells_gex_scREA_rep_A_B_2B_20220610.txt",
            sep="\t",row.names=FALSE,quote=FALSE)


