# getzlab-pvalue_combination
p-value combination for PCAWG drivers

Note: we are using parts of the original code for the Empirical Brown's method (Poole et al. Bioinformatics, 2016). For details, check the following github repo: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM

The version of R is 3.4.

Please install the following package before running the command: reshape

The inputs to the R script are: tissue, target, dir, sif_filepath.
The outputs contain a a diagnostic QQ plot of observed vs. expected p-values for all input p-value sets, a table with combined p-values, a report that indicates which p-value sets were included and a table of combined p-values after removing some outlier methods (.combined_p_values.automatic_method_removal.txt). For details please see XXX

Steps:

1. Copy ebm.R from https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM into pvalue_combination/

2. Unarchive example_data/PCAWG_meta_Carcinoma_cds_pvalues.zip

3. Run the following command to integrate p values for Carcinoma CDS (protein-coding) regions:

Rscript pvalue_combination/run_pvalue_combination_comparison_automatic_removal.R meta_Carcinoma CDS output/ example_data/SIF_input.txt

The outputs will be saved under output/ 

Some explanations for columns in .automatic_method_removal_report.txt:

-median_sig:

-upper_thresh:

-removed_for_inflation:
