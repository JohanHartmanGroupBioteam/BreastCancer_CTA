## extract scaled expr obj

#### loading PKGs ####
library(ggplot2)
library(tidyverse)

library(futile.logger)

library(infercnv)
#library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)


library(SpatialInferCNV)

library(Seurat)
library(dplyr)
library(magrittr)

library(dendextend)

#### scale cnv  ####

scale_infercnv = function(infercnv_obj, x.range = "auto"){
  
  plot_data = infercnv_obj@expr.data
  x.center= mean(plot_data)
  if (! any(is.na(x.range))) {
    
    
    if ( (length(x.range) == 1) & (x.range[1] == "auto") ) {
      
      # examine distribution of data that's off-center, since much of the center could
      # correspond to a mass of data that has been wiped out during noise reduction
      quantiles = quantile(plot_data[plot_data != x.center], c(0.01, 0.99))
      
      # determine max distance from the center.
      delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
      low_threshold = x.center - delta
      high_threshold = x.center + delta
      x.range = c(low_threshold, high_threshold)
      
      flog.info(sprintf("plot_cnv(): auto thresholding at: (%f , %f)", low_threshold, high_threshold))
      
    } else {
      
      # use defined values
      low_threshold = x.range[1]
      high_threshold = x.range[2]
      
      if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
        stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
      }
    }
    
    plot_data[plot_data < low_threshold] <- low_threshold
    plot_data[plot_data > high_threshold] <- high_threshold
    
    infercnv_obj@expr.data <- plot_data  #because used again below...
    
  }
  
  return(infercnv_obj)
}


#### loading data ####

# respath="/Users/qiaoyang/MountScilife/projects/3.CIIR/2.benchmark/run/inferCNV/res_spinfercnv"

respath="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/inferCNV/res_spinfercnv"

files = list.files(path=respath,pattern = "run.final.infercnv_obj"  ,full.names = TRUE, recursive = TRUE )

for (file in files){
  
  print(file)
  
  infercnv_obj = readRDS( file )
  
  infercnv_obj = scale_infercnv(infercnv_obj) ## do scale and clean data
  
  saveRDS(infercnv_obj, file =  str_replace(file, "run.final.infercnv_obj","run.final.scaled.infercnv_obj" ) )
  
}


