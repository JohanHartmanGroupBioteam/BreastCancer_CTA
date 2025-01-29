# BreastCancer_CTA
R codes used for analysis of results from computational tissue annotation (CTA) of breast cancer sample

- **ST_deconvolution.Rmd**: code used to prepare single-cell RNAseq reference matrix and spatial transcriptomics matrix for deconvolution

- **geojson_tissue_level_correlation.ipynb**: code used to align visium spots and tissue-level CTA results. After exporting the tissue-level CTA results as geojson files and creating visium spots using code in ST_deconvolution.Rmd, this code is applied to calculate the overlay percentage for tumor and stromal+immune compartments 

- **Digital_pathology_deconvolution_correlation.Rmd**: code used to perform Spearman's correlation between CTA results and deconvolution results at both cell- and tissue-level

- **ST_analysis.Rmd**: code used for spatial transcriptomic data integration, analysis, and visualization