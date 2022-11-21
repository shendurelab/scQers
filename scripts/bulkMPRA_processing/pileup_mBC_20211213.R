

library(stringr)
require(dplyr)

args <- commandArgs(trailingOnly=TRUE)
in_file_oi <- args[1]
out_file_oi <- args[2]

print("Reading table")
df <- read.table(in_file_oi,sep="\t",header=TRUE)

print("Processing and piling up reads")
BC_UMI <- df$BC
split_BC_UMI <- str_split(BC_UMI,"_")
mBC <- unlist(lapply(split_BC_UMI,"[[",1))
UMI <- unlist(lapply(split_BC_UMI,"[[",2))
df$mBC <- mBC
df$UMI <- UMI
df_noG <- df %>% filter(! (UMI=="GGGGGGGGGG") )
df_noG$unique_mBC_UMI <- rep(1,dim(df_noG)[1])

quant_mBC <- df_noG %>% group_by(mBC) %>% summarise(n_UMI_per_mBC = sum(unique_mBC_UMI), 
                                                    n_reads_per_mBC = sum(BC_counts))

print("Writing output")
write.table(quant_mBC,out_file_oi,
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
