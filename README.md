# scQers 🍢
Scripts, plasmid maps, and amplicon maps for: Multiplex profiling of developmental enhancers with quantitative, single-cell expression reporters

Note: scripts contain herein do not currently constitute streamlined pipelines, and are shared for transparency. Iterations and improvement will be made in near future. 

Organization of repositories: 

Plasmid maps can be found in "plasmid_maps" (refer to manuscript Methods for further details). 

Structure of custom sequencing amplicons used in this study are listed in "custom_amplicon_structures". 

Scripts is separated in subdirectories organized by categories as detailed below with short descriptive of contents. 


## scQer_analysis

### pre_processing
Gene expression count matrices (from output of cellranger) are filtered (gene with more than 3 cells detected, cell barcodes with at least 50 genes detected, total GEx UMI counts above threshold [sample/depth specific, but typically 500-1000 UMI/cell barcode], and mitochondrial fraction within 1 to 15%). The filtered cell barcodes have their scrublet score determined (GEx_compute_scrublet_score.py). After filtering for cell barcodes with scrublet score below the bimodal threshold (determined from simulated doublets), replicates are merged, dimensionally reduced, and clustered (e.g., GEx_scQer_processing.R). Sub-clustering of identified clusters is performed to further remove likely doublet cell barcodes (subclustering_GEx_cluster_doublet_cleanup.R). oBC and mBC fastqs are converted to raw cell barcode and BC UMI count (error corrected) table (respectively: oBC_scQer_preprocessing.sh and mBC_scQer_preprocessing.sh, see example outputs: mBC_A1_get_bc_v3_no_G_cleaned_UMI_20220525.txt, oBC_mEB_A1_get_bc_v3_no_G_cleaned_UMI_20220604.txt). These were then filtered for valid cell barcodes (as determined from the gene expression data) and valid BC (as predetermined from the oBC-CRE-mBC subassembly of the scQer library). The filtered data was then be merged (see GSE217686: 
GSE217686_assigned_oBC_CRE_mBC_joined_counts_sc_rep_mEB_series.txt.gz, 
GSE217689: GSE217689_assigned_oBC_CRE_mBC_joined_counts_sc_rep_promoter_series.txt.gz). 

### statistical_testing
Functions bootstrap_scQer_activity_quantification.R and permutation_scQer_specificity_quantification.R are used to generate the bootstrap (activity) and permutation (specificity) resamplings (these are run through arrayed jobs, generating a file for each CRE and replicate, which are then combined to a single file with  combine_bootstraps.R). The empirical corrected p-values are generated in compute_empirical_p_values.R, which are used to identify reproducibly and significantly active/specific CREs following criteria listed in Methods. 

### clonotype_analysis
oBC_clonotype_assignment_fisher_exact.R takes the cell barcode/oBC UMI count table, and groups cells based on co-detection by Fisher exact test (method inspired by Wang et al, BMC Genomics, 2022), generating a raw clonotype set. Script refine_clonotype.R then filters through the raw groups of cells, leading to a high confidence set of clonotypes. In assign_cell_to_clonotype.R, these refined clonotypes (set of co-detected oBC) are then used to assign cells back  with low stringency of the fraction of detected oBC (>50%), but a high stringency to not detect barcodes from other clones (to filter out doublets), generate_final_clonotype_lists.R collates the information in final metadata tables. oBC_detection_precision_recall.R performs the analysis of oBC capture based on the groud-truth clones identified.


## plasmid_library_subassembly:

### oBC-mBC subassembly:
Preprocessing is performed with (arr_mBC_oBC_association_preprocessing.sh). Given that the result is to count barcodes, preprocessing for oBC-mBC association uses a number of the same dependency as bulk MPRA preprocessing. The raw count table is then processed by thresholding and determining likely non-uniquely paired oBC-mBC in mBC_oBC_processing.R. Example final output can be downloaded from GEO (GSE217681, GSE217681_p025_recloned_complex_mBC_oBC_subassembly.txt.gz)

### oBC-CRE subassembly:
Preprocessing script (oBC_CRE_association_preprocessing.sh) is run to generate a barcode to CRE count table, with a tally of read counts, orientation, and tagmentation positions (see example_out: p55_1x_10m2_oBC_devCRE_suba_S2_full_pileup_k2_20220818.txt.gz). Script generate_oBC_CRE_association.R then filters further for bona fide associations based on expected position of the aligned reads, unicity of oBC to CRE pairing, and read counts (example_outs: p55_1x_10m2_oBC_devCRE_suba_S2_w_multimapper_final_triplets_20220818.txt). This table is further filtered to the final subassembly by combining with the promoter series triplet and retaining unique mBC (see GEO: GSE217681, GSE217681_final_subassembly_oBC_devCRE_mBC_triplets_mEB_series.txt.gz).


## bulkMPRA_processing: 
The bash pipeline (arr_bulkMPRA_preprocessing.sh) is first run, generated a raw count table (see BC_quant file in example_output directory). Then the counts from RNA and DNA data for each barcodes are associated with pre-determined matched regulatory elements (from subassembly stage) using script generate_count_table_w_CRE_bulkMPRA.R (see example_output: bulk_MPRA_mEB_batch1_seq027_v2_20220826.txt).


## data_integration
Minimal Seurat objects from the scRNA-seq data were created (from scQer day 21 mEB data,  in vivo Pijuan-Sala et al, Nature, 2019), and integrated following Butler et al, Cell 2019. Following integration, the cell-type labels from the 100 closest neighbours from Pijuan-Sala to each cell from the mEB dataset were saved (combine_integrate_scRNA_obj.R). The top 10 nearest neighbors’ labels were transfered if >=60% were of the same type (otherwise uncertain) using script aggregate_transfered_labels.R. 

## TFBS_analysis
TFBS affinity for 8-mers were downloaded from Uniprobe (http://the_brain.bwh.harvard.edu/uniprobe/) for endodermal transcription factors Gata4, Foxa2, and Sox17 and converted to relative affinities (essentially: [value-mode]/[max-mode], see table: TFBS_relative_affinity_tables_endo_TFs_Uniprobe.txt). Sequences around loci considered were then scanned for TFBS above a given threshold using TFBS_locus_scan.R (example output: TFBS_hits_Bend5_100kb_pm_TSS.txt). The results were then aggregated to 500 bp windows in distribution_TFBS_counts.R and compared to the TFBS counts found in distal CREs. 



