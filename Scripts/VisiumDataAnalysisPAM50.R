#' ---
#' title: "Visium BCSAx data analysis & PAM50 subtyping (naive spatial minibulk)"
#' author: "Emmanouil G. Sifakis"
#' license: "CC BY 4.0"
#' date: "`r format(Sys.Date(), format = '%B %d, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      theme: cosmo
#'      toc_float: true
#'      number_sections: true
#' ---
#' ```{r global_options, include=FALSE}
#' knitr::opts_chunk$set(
#' fig.width=10,
#' fig.height=5,
#' fig.path='figures/',
#' dev='svg',
#' echo=F,
#' warning=FALSE,
#' message=FALSE,
#' knitr.kable.NA = '')
#' ```

rm(list=ls()) # clear all

#' # Load packages
#' 
library(tidyverse)
library(Seurat)
suppressPackageStartupMessages(library(STutility))


#' # Set paths
#' 
data_dir <- "P:/Hartman group/CIIR"
output_dir <- "C:/Users/emmsif/OneDrive - Karolinska Institutet/OnkPat/projects/JohanHartman/CIIR"


#' # Load filtered ST data
#' 

files <- list.files(data_dir, pattern = "filtered_feature_bc_matrix.h5",recursive = T, full.names = T)


sample_info <- readxl::read_xlsx(file.path(output_dir,"ST-AR_Visium_samples_se_QuPathValidation_updated.xlsx"),sheet = "Sheet1")
colnames(sample_info) <- c("pid","tid","vid")

for (f in 1:length(files)) {

  fid <- files[f]
  
  
  sid <- sample_info[sample_info$vid %in% sub("/outs/filtered_feature_bc_matrix.h5","",sub(paste0(data_dir,"/"),"",fid)),]
  vid <- sid$vid
  sid <- paste0(sid$pid,sid$tid)
  
  
  load(file=file.path(output_dir,paste0(sid,".RData")))
  
  
  #' # Analyzing sample
  #' 
  knitr::kable(vid)
  knitr::kable(sid)
  
  #' # Analysis
  #' 
  
  #' ## Normalization
  #' 
  
  bcsa <- SCTransform(bcsa, assay = "Spatial", verbose = F, method = "poisson")
  
  #' ## Plot gene expression of individual genes
  #' 
  #+ fig.width=10, fig.height=30
  SpatialFeaturePlot(bcsa, features = c("ESR1", "PGR", "ERBB2", "MKI67"))
  
  #' # PAM50 subtyping (naive spatial minibulk)
  #'
  
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(genefu)
  library(RColorBrewer)
  
  data("pam50", package = "genefu")
  rownames(pam50$centroids.map) <- pam50$centroids.map$EntrezGene.ID
  
  lkup <- c(
    "LumA" = "LA",
    "LumB" = "LB",
    "Her2" = "H2",
    "Basal" = "BL",
    "Normal" = "NBL"
  )
  
  subtype_lkup <- c(
    "LA" = "Luminal A",
    "LB" = "Luminal B",
    "H2" = "HER2-enriched",
    "BL" = "Basal-like",
    "NBL" = "Normal breast-like")
  
  colsubtypecd <- c(
    "LA" = "#2a3188",
    "LB" = "#419ad2",
    "H2" = "#d4279c",
    "BL" = "#97191e",
    "NBL" = "#66c530",
    "NaN" = "black"
  )

  
  #' ## Add ENTREZID as meta data
  #'
  
  # To the Spatial assay
  
  genes <- rownames(bcsa@assays$Spatial) %>% as.data.frame()
  
  colnames(genes) <- "Gene.Symbol"
  genes$EntrezGene.ID <- mapIds(org.Hs.eg.db,
                                keys=genes$Gene.Symbol,
                                column="ENTREZID",
                                keytype="SYMBOL",
                                multiVals="first")
  
  genes[which(genes$Gene.Symbol=="C5orf30"),]$EntrezGene.ID <- "90355"
  genes[which(genes$Gene.Symbol=="MKL2"),]$EntrezGene.ID <- "57496"
  genes[which(genes$Gene.Symbol=="GARS"),]$EntrezGene.ID <- "2617"
  genes[which(genes$Gene.Symbol=="FAM49B"),]$EntrezGene.ID <- "51571"

  rownames(genes) <- genes$Gene.Symbol
  
  
  bcsa@assays$Spatial@meta.data$Gene.Symbol <- genes$Gene.Symbol
  bcsa@assays$Spatial@meta.data$EntrezGene.ID <- genes$EntrezGene.ID
  
  #' ## (naive spatial minibulk) AIMS: Absolute Assignment of Breast Cancer Intrinsic Molecular Subtype
  #'
  #' Notes: AIMS [Paquet and Hallett, 2015](https://academic.oup.com/jnci/article/107/1/dju357/900920), on the raw counts
  
  AIMS <- molecular.subtyping(sbt.model="AIMS",
                              data=t(bcsa@assays$Spatial$counts) %>% as.data.frame(),
                              annot=bcsa@assays$Spatial@meta.data,
                              do.mapping=T,
                              verbose = T)
  
  subtypecd_aims <- factor(AIMS$subtype, levels = c("LumA","LumB","Her2","Basal","Normal"))
  names(subtypecd_aims) <- rownames(AIMS$subtype)
  subtypecd_aims <- plyr::revalue(subtypecd_aims, lkup)
  
  ### Append results
  subtypecd_aims <- subtypecd_aims[match(Cells(bcsa), names(subtypecd_aims))]
  bcsa$subtypecd_aims <- subtypecd_aims
  
  
  #' ### Plot (naive spatial) AIMS-PAM50 onto the tissue section
  #'
  
  SpatialDimPlot(bcsa, group.by = c("subtypecd_aims"), cols=colsubtypecd)

  
  #' ### Export results
  #'
  save(bcsa, file=file.path(output_dir,paste0(sid,"_processed.RData")))
  rm(bcsa)
  gc()

}
