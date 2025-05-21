# BreastCancer_CTA
R codes used for analysis of results from computational tissue annotation (CTA) of breast cancer sample

- **VisiumDataQC.R**: code used to perform quality control of the Visium data (including mitochondrial, ribosomal, and hemoglobin percent calculation)

- **ST_deconvolution.Rmd**: code used to prepare single-cell RNAseq reference matrix and spatial transcriptomics matrix for deconvolution

- **Xenium_CTA_comparison.Rmd**: code used to compare output from CTA and cell type markers' expression from Xenium platform. 

- **Digital_pathology_deconvolution_correlation.Rmd**: code used to perform Spearman's correlation between CTA results and deconvolution results at both cell- and tissue-level

- **VisumDataAnalysisPAM50.R**: code used for performing PAM50 analysis using AIMS method on spatial transcriptomic data

- **ST_analysis.Rmd**: code used for spatial transcriptomic data integration, analysis, and visualization

- **SpatialVDJ_visualization.Rmd**: code used for visualization of immune clones identified using SpatialVDJ

