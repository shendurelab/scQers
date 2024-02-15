library(tidyverse)
library(ggplot2)

# function argument
args <- commandArgs(trailingOnly=TRUE)
file_BC_in <- args[1]
file_BC_out <- args[2]
read_thresh <- args[3] %>% as.numeric
fig_out_name <- args[4]

# load BC table
df_BC <- read.table(file_BC_in,header=TRUE)

# count reads per BC
df_BC_count <- df_BC %>% group_by(BC) %>% summarize(BC_read_counts=length(BC))

# select BC in the high count mode
df_BC_hi <- df_BC_count %>% filter(BC_read_counts>read_thresh)

# print list of high count BCs
write.table(df_BC_hi,file_BC_out,
            row.names=FALSE,sep="\t",quote=FALSE)

# plot distribution of counts
plt_pdf_BC_counts <- ggplot(df_BC_count) + stat_bin(aes(x=BC_read_counts),geom="step",bins=50)+
  geom_vline(aes(xintercept=read_thresh),linetype="dashed",color="orange")+
  scale_x_log10()+scale_y_log10()+
  labs(x="Read counts per BC",y="Distribution (number of BC)")

pdf(fig_out_name,width=3,height=2.5)
print(plt_pdf_BC_counts)
dev.off()

