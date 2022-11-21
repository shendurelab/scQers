import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

file_names = ["gex_2B1_filtered_10x_mtx_20220610","gex_2B2_filtered_10x_mtx_20220610"]
starting_dir = "/Users/jbl/Documents/UW/Data/sequencing_runs/seq031_bottlenecked_mEB_scREA_20220607/GEX_outs_20220608/"

for file_name in file_names:
    
    print("processing "+file_name)
    dir_oi = starting_dir+file_name
    
    cell_bc = pd.read_csv(dir_oi+"/barcodes.tsv",header=None)
    counts_matrix = scipy.io.mmread(dir_oi+"/matrix.mtx")
    genes = np.array(scr.load_genes(dir_oi+"/genes.tsv", column=1))

    print(len(genes))
    print(counts_matrix.shape)

    scrub_z = scr.Scrublet(np.transpose(counts_matrix))
    doublet_scores_z, predicted_doublets_z = scrub_z.scrub_doublets(n_prin_comps=30, 
                                                                mean_center=True, 
                                                                normalize_variance=True)
    df=pd.DataFrame()
    df['cell_bc']=cell_bc[0]
    print(len(cell_bc[0]))
    print(len(doublet_scores_z))
    df['doublet_scores_z'] = doublet_scores_z
    
    df.to_csv(dir_oi+'/scrublet_score_scRNA'+file_name+'.txt',sep='\t',index=False)                                                               
                                                                
    #print(le))
                                                                
    fig, axs = scrub_z.plot_histogram()
    print(fig)
    fig.savefig(dir_oi+"/doublet_score_histograms_scrublet_" + file_name + ".pdf")
