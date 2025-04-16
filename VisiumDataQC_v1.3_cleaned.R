#' ---
#' title: "Visium BCSAx data Quality Control"
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


#' # Set paths
#' 
data_dir <- "P:/Hartman group/CIIR"
output_dir <- "C:/Users/emmsif/OneDrive - Karolinska Institutet/OnkPat/projects/JohanHartman/CIIR"


#' # Load ST data (and filter for spots that does not overlap with the tissue)
#' 

files <- list.files(data_dir, pattern = "filtered_feature_bc_matrix.h5",recursive = T, full.names = T)


sample_info <- readxl::read_xlsx(file.path(output_dir,"ST-AR_Visium_samples_se_QuPathValidation_updated.xlsx"),sheet = "Sheet1")
colnames(sample_info) <- c("pid","tid","vid")
  

for (f in 1:length(files)) {

  fid <- files[f]
  
  sid <- sample_info[sample_info$vid %in% sub("/outs/filtered_feature_bc_matrix.h5","",sub(paste0(data_dir,"/"),"",fid)),]
  vid <- sid$vid
  sid <- paste0(sid$pid,sid$tid)
  
  bcsa <- Load10X_Spatial(
    data.dir = dirname(fid),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sid,
    filter.matrix = TRUE
  )
  bcsa$orig.ident <- sid
  bcsa$type <- "tumor"
  
  
  
  
  ## Crop folded tissue in region "BCSA1TumA2"
  
  if (sid=="BCSA1TumA2"){
  
    # Visualize the original region
      Seurat::SpatialDimPlot(bcsa, crop = F,  alpha = c(0.5, 1))
    
    # coords <- GetTissueCoordinates(bcsa)
    coords <- GetTissueCoordinates(bcsa, scale = "hires", cols = c("imagerow", "imagecol"))
    coords$x %>% summary
    coords$y %>% summary
    
    # Coordinates for filtering cells that belong to the folded tissue area
  
    x_min <- 344.8
    x_max <- 1110
    y_min <- 900
    y_max <- 1746.8 
    
    # Cells to remove that belong to the folded tissue area
    cells_to_remove <- rownames(coords[coords$x >= x_min & coords$x <= x_max & coords$y >= y_min & coords$y <= y_max, ])
    
    # Subset the Seurat object based on the filtered coordinates
    filtered_bcsa <- subset(bcsa, cells = setdiff(Cells(bcsa), cells_to_remove), invert=F)
    
    # Visualize the cropped region
    Seurat::SpatialDimPlot(filtered_bcsa, crop = F,  alpha = c(0.5, 1))
  
    bcsa <- filtered_bcsa
  }
  
  
  
  #' ## Import Space Ranger graph-based clustering (https://www.10xgenomics.com/support/software/space-ranger/algorithms-overview/gene-expression)
  #' 
  
  clusters <- readr::read_csv(file.path(dirname(fid),"analysis/clustering/graphclust/clusters.csv")) %>% as.data.frame()
  rownames(clusters) <- clusters$Barcode
  
  cluster_numbers <- clusters[Cells(bcsa),]$Cluster
  names(cluster_numbers) <- colnames(bcsa)
  bcsa <- AddMetaData(
    object = bcsa,
    metadata = cluster_numbers,
    col.name = 'spaceranger_clusters'
  )
  # head(bcsa[[]])
  
  
  
  #' # Analyzing sample
  #' 
  
  knitr::kable(vid)
  knitr::kable(sid)
  
  
  #' # Quality control
  #' 
  
  bcsa <- PercentageFeatureSet(bcsa, "^MT-", col.name = "percent_mito")
  bcsa <- PercentageFeatureSet(bcsa, "^RP[SL]", col.name = "percent_ribo")
  bcsa <- PercentageFeatureSet(bcsa, "^HB[^(P)]", col.name = "percent_hb")
  
  
  feats <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribo", "percent_hb")
  
  #+ fig.width=10, fig.height=10
  VlnPlot(bcsa, split.by = "type", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
  
  #+ fig.width=10, fig.height=30
  SpatialFeaturePlot(bcsa, features = feats)
  
  
  
  
  #' ## Filter
  #' 
  
  bcsa <- bcsa[, bcsa$percent_mito < 25 & bcsa$percent_hb < 20 & bcsa$nFeature_Spatial > 500]
  
  #+ fig.width=10, fig.height=30
  SpatialFeaturePlot(bcsa, features = feats)
  
  
  #' ## Top expressed genes
  #' 
  
  C = bcsa@assays$Spatial$counts
  C@x = C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
  #+ fig.width=10, fig.height=10
  boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.5, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
  
  
  #' ## Filter genes
  #' 
  
  dim(bcsa)
  
  # Filter Mitochondrial genes
  bcsa <- bcsa[!grepl("^MT-", rownames(bcsa)), ]
  
  dim(bcsa)
  
  
  #' ## Top expressed genes (after filtering not interesting genes)
  #' 
  
  C = bcsa@assays$Spatial$counts
  C@x = C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
  #+ fig.width=10, fig.height=10
  boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.5, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
  
  
  #' ### Export results
  #'
  save(bcsa, file=file.path(output_dir,paste0(sid,".RData")))
  rm(bcsa)
  gc()


}
