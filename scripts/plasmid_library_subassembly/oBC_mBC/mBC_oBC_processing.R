require(dplyr); require(data.table); require(ggplot2);
require(cowplot); theme_set(theme_cowplot()); require(stringr)
library(scales)
library(Biostrings)
library(patchwork)
library(stringdist)

library(tidyverse)
library(tictoc)

# start processing the list of barcodes

dir_out <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq020_p25_asso_p22_mBC_sc_rep_CS2_v2_bulkoBC_sci_enrich/p25_mBC_oBC"
setwd(dir_out)


# cdf_BC_counts <- ecdf(df_oBC_mBC_2$oBC_mBC_read_counts)
# 
# xplot <- 10^(seq(-0.3,6,0.001))
# xlabels <- c(1,10,100,1000,10000,1E5)
# xlims <- c(0.5,3E5)
# ylims <- c(1E-5,1)
# ylabels <- c(1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1)
# h_fig <- ggplot() + geom_step(aes(x=xplot,y=1-cdf_BC_counts(xplot)),color='black')+
#     scale_x_log10(breaks = xlabels,
#                   labels = trans_format("log10", math_format(10^.x)),
#                   limits = xlims)+
#     scale_y_log10(breaks = ylabels,
#                   labels = trans_format("log10", math_format(10^.x)),
#                   limits = ylims)+annotation_logticks()+
#     geom_vline(xintercept=count_thresh[ind],color="red",linetype="dashed")+coord_fixed()+
#     labs(x='Read counts [unique mBC-oBC pairs]',y='1-CDF',title=file_names[ind])
# h_fig
# 
# sum(df_oBC_mBC_2$oBC_mBC_read_counts>count_thresh[ind])


## axis limits for QC plots
# xlims_prop <- c(0.001,2)
# xlabels_prop <- c(1E-3,1E-2,1E-1,1)
# 
# ylims_prop <-  c(0.001,2)
# ylabels_prop <- c(1E-3,1E-2,1E-1,1)


dir_data <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq020_p25_asso_p22_mBC_sc_rep_CS2_v2_bulkoBC_sci_enrich/p25_mBC_oBC/condensed/"
file_names <- c("p25_mBC_oBC_asso_0p001x_S17", "p25_mBC_oBC_asso_0p01x_S16","p25_mBC_oBC_asso_0p1x_S15", "p25_mBC_oBC_asso_1x_S14")
suffix <- "_pear_condensed_20220320.txt.gz"
count_thresh <- c(300,30,12,4) # determined by inspection of the distribution of barcode couts 
prop_thresh <- 0.95

for (ind in seq(length(file_names))){
    

    file_oi <- paste0(dir_data,file_names[ind],suffix)
    
    df <- read.table(file_oi, sep="\t", header=TRUE)
    
    oBC_mBC <- df$BC
    split_oBC_mBC <- str_split(oBC_mBC,"_")
    oBC <- unlist(lapply(split_oBC_mBC,"[[",1))
    mBC <- unlist(lapply(split_oBC_mBC,"[[",2))
    
    df_oBC_mBC <- data.frame(oBC,mBC,oBC_mBC_read_counts=df$BC_counts) 
    
    bool_mBC_N <- str_detect(df_oBC_mBC$mBC,'N')
    bool_oBC_N <- str_detect(df_oBC_mBC$oBC,'N')
    df_mBC_oBC_no_N <- df_oBC_mBC %>% filter(!bool_mBC_N & !bool_oBC_N)
    df_mBC_oBC_no_N_no_Gs <- df_mBC_oBC_no_N %>% filter(mBC!="GGGGGGGGGGGGGGG" & oBC!="GGGGGGGGGGGGGGGG")
    
    df_oBC_mBC_2 <- df_mBC_oBC_no_N_no_Gs
    
    ## QC plots and  distribution
    # cdf_BC_counts <- ecdf(df_oBC_mBC_2$oBC_mBC_read_counts)
    # 
    # xplot <- 10^(seq(-0.3,6,0.001))
    # xlabels <- c(1,10,100,1000,10000,1E5)
    # xlims <- c(0.5,3E5)
    # ylims <- c(1E-5,1)
    # ylabels <- c(1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1)
    # h_fig <- ggplot() + geom_step(aes(x=xplot,y=1-cdf_BC_counts(xplot)),color='black')+
    #     scale_x_log10(breaks = xlabels,
    #                   labels = trans_format("log10", math_format(10^.x)),
    #                   limits = xlims)+
    #     scale_y_log10(breaks = ylabels,
    #                   labels = trans_format("log10", math_format(10^.x)),
    #                   limits = ylims)+annotation_logticks()+
    #     geom_vline(xintercept=count_thresh[ind],color="red",linetype="dashed")+
    #     labs(x='Read counts [unique mBC-oBC pairs]',y='1-CDF',title=file_names[ind])
    
    # h_fig2 <- ggplot() + geom_step(aes(x=xplot,y=cdf_BC_counts(xplot)),color='black')+
    #     scale_x_log10(breaks = xlabels,
    #                   labels = trans_format("log10", math_format(10^.x)),
    #                   limits = xlims)+ylim(c(0.3,1))+
    #     geom_vline(xintercept=count_thresh[ind],color="red",linetype="dashed")+
    #     labs(x='Read counts [unique mBC-oBC pairs]',y='CDF')
    
    # fig_name <- sprintf("p25_v2_distribution_oBC_mBC_counts_1mCDF_%s_20220320.pdf",file_names[ind])
    # pdf(fig_name, width=4, height=4)
    # print(h_fig)
    # dev.off()
    # 
    # fig_name <- sprintf("p25_v2_distribution_oBC_mBC_counts_CDF_%s_20220320.pdf",file_names[ind])
    # pdf(fig_name, width=4, height=4)
    # print(h_fig2)
    # dev.off()
    
    
    
    df_hi_mBC_oBC_pairs <- df_oBC_mBC_2 %>% filter(oBC_mBC_read_counts>count_thresh[ind])
    df_hi_mBC_oBC_pairs <- df_hi_mBC_oBC_pairs %>% transform(prop_mBC_oBC_counts = oBC_mBC_read_counts/sum(oBC_mBC_read_counts))
    
    
    # calculating proportion of reads to a given mBC or oBC. 
    df_sum_mBC <- df_hi_mBC_oBC_pairs %>% group_by(mBC) %>% summarise(sum_counts_per_mBC = sum(oBC_mBC_read_counts))
    df_hi_mBC_oBC_pairs2 <- left_join(df_hi_mBC_oBC_pairs,df_sum_mBC)
    df_sum_oBC <- df_hi_mBC_oBC_pairs2 %>% group_by(oBC) %>% summarise(sum_counts_per_oBC = sum(oBC_mBC_read_counts))
    df_hi_mBC_oBC_pairs3 <- left_join(df_hi_mBC_oBC_pairs2,df_sum_oBC)
    
    df_hi_mBC_oBC_pairs <- df_hi_mBC_oBC_pairs3 %>% transform(prop_oBC = oBC_mBC_read_counts/sum_counts_per_oBC, prop_mBC = oBC_mBC_read_counts/sum_counts_per_mBC)
    
    # # QC plots 
    # print(file_names[ind])
    # print(sprintf("# pairs with high read counts (>%d):",count_thresh[ind]))
    # print(dim(df_hi_mBC_oBC_pairs)[1])
    # 
    # print("prop unique pairs:")
    # print(mean(df_hi_mBC_oBC_pairs$prop_oBC>prop_thresh & df_hi_mBC_oBC_pairs$prop_mBC>prop_thresh))
    # 
    # prop_fig1 <- ggplot(df_hi_mBC_oBC_pairs) + geom_step(aes(x=prop_oBC),stat="ecdf")+theme_gray() + geom_vline(xintercept=prop_thresh,color="red",linetype="dashed")+
    #     labs(y="CDF")+coord_fixed()
    # prop_fig2 <- ggplot(df_hi_mBC_oBC_pairs) + geom_step(aes(x=prop_mBC),stat="ecdf")+theme_gray() + geom_vline(xintercept=prop_thresh,color="red",linetype="dashed")+
    #     labs(y="CDF")+coord_fixed()
    # prop_fig3 <- ggplot(df_hi_mBC_oBC_pairs) + geom_step(aes(x=prop_oBC),stat="ecdf")+theme_gray() +xlim(c(0.9,1))+ geom_vline(xintercept=prop_thresh,color="red",linetype="dashed")+
    #     labs(y="CDF")
    # prop_fig4 <- ggplot(df_hi_mBC_oBC_pairs) + geom_step(aes(x=prop_mBC),stat="ecdf")+theme_gray() +xlim(c(0.9,1))+ geom_vline(xintercept=prop_thresh,color="red",linetype="dashed")+
    #     labs(y="CDF")
    # 
    # 
    # h_fig3 <- ggplot(df_hi_mBC_oBC_pairs)+geom_bin_2d(aes(x=prop_oBC,y=prop_mBC),n=100)+
    #     scale_x_log10(breaks = xlabels_prop,
    #                   labels = trans_format("log10", math_format(10^.x)),
    #                   limits = xlims_prop)+
    #     scale_y_log10(breaks = ylabels_prop,
    #                   labels = trans_format("log10", math_format(10^.x)),
    #                   limits = ylims_prop)+annotation_logticks()+coord_fixed()+
    #     labs(x='Proportion reads to oBC',y='Proportion reads to mBC')
    # fig_name <- sprintf("p25_v2_prop_reads_oBC_mBC_%s_20220320.pdf",file_names[ind])
    # pdf(fig_name, width=5.5, height=4)
    # print(h_fig3)
    
    # plt_design <- "AAB
    #                CD#"
    # 
    # fig_name <- sprintf("p25_v2_oBC_mBC_asso_summary_%s_20220320.pdf",file_names[ind])
    # pdf(fig_name, width=11, height=12)
    # 
    # fig_final <- h_fig+h_fig3+prop_fig1+prop_fig2+plot_layout(design=plt_design)
    # print(fig_final)
    # dev.off()
    # 
    
    df_hi_mBC_oBC_pairs_2 <- df_hi_mBC_oBC_pairs %>% transform(bool_prop_mBC=prop_mBC>prop_thresh, bool_prop_oBC=prop_oBC>prop_thresh, bool_valid=(prop_oBC>prop_thresh & prop_mBC>prop_thresh))
    
    # write output file: 
    file_name_table <- sprintf("final_table_p25_v2_oBC_mBC_subassembly_%s_20220320.txt",file_names[ind])
    write.table(df_hi_mBC_oBC_pairs_2, file_name_table,
                quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
    
    
}
   

