# getzlab-pvalue_combination
p value combination for PCAWG drivers

Note: we are using parts of the original code for Empirical Brown's method (Poole et al. Bioinformatics, 2016). For details, check the following github repo: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM

The version of R is 3.4.

The inputs to the R script are: tissue, target, dir, sif_filepath.
The outputs contain a p values QQ plot, a table with combined p values, a report with methods to be removed and an integrated table of combined p values after removing some outlier methods.

Steps:

1. Copy ebm.R from https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM into pvalue_combination/

2. Unarchive input_pvalues/PCAWG_meta_Carcinoma_cds_pvalues.zip

3. Run the following command to integrate p values for Carcinoma CDS region:

Rscript pvalue_combination/run_pvalue_combination_comparison_automatic_removal.R meta_Carcinoma CDS output/ example_data/SIF_input.txt

The outputs will be saved under output/
