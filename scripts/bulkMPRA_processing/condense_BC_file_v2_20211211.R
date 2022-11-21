

args <- commandArgs(trailingOnly=TRUE)
in_file_oi <- args[1]
out_file_oi <- args[2]
var1_to_pile <- args[3]
var2_to_pile <- args[4]

print("Reading table")
df_BC <- read.table(in_file_oi, sep = '\t', header=TRUE)

print("combining two components")
BC_UMI <- paste0(df_BC[[var1_to_pile]], "_", df_BC[[var2_to_pile]])


print("Piling up reads")
BC_UMI_counts <- table(BC_UMI)

print("Creating pile-up df")
df_BC_UMI_counts <- data.frame(BC=names(BC_UMI_counts),
                                BC_counts=as.numeric(BC_UMI_counts))

print("Writing output")
write.table(df_BC_UMI_counts,out_file_oi,
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')