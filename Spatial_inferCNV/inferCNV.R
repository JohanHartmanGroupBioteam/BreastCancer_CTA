#### loading libs ####

## run in r-environment 
## conda activate r-environment

library(zellkonverter)
library(anndata)
library(ggplot2)
library(tidyverse)

library(infercnv)
#library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)

library(SpatialInferCNV)

library(zellkonverter)
library(anndata)
library(ggplot2)
library(tidyverse)

library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)

library(SpatialInferCNV)

library(magrittr)

library(dplyr)
library(magrittr)

#### loading expr data #### 
path="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/inferCNV/merged_data"
respath="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/inferCNV/res_spinfercnv"

#path="/Users/qiaoyang/MountScilife/projects/3.CIIR/2.benchmark/run/inferCNV/merged_data"

print("Load BCSA1")
load(paste( path, "BCSA1_merge.Rdata", sep = "/"))

# Layers:
#   counts.1.1.1, counts.2.1.2, counts.1.2.3, counts.2.2.4
# BCSA1TumA1, BCSA1TumA2, BCSA1TumB1, BCSA1TumB2


# filtered 
# bcsa = subset(bcsa1, subset = nFeature_Spatial > 500 & percent_mito < 25 & percent_hb < 0.2)

# counts = bcsa1@assays$Spatial$counts.1.1
# counts = SeuratObject::LayerData(object =bcsa1[["Spatial"]], layer = "counts.1.1" )


#### load gene annotation ####

# library(rtracklayer)
## on local for convience
# GTFfile="/Users/qiaoyang/MountScilife/projects/3.CIIR/2.benchmark/run/inferCNV/merged_data/gencode.v32.annotation.gtf.gz"
# GTF = rtracklayer::import(GTFfile)
# GTF = data.frame(GTF)
# save(GTF, file = "/Users/qiaoyang/MountScilife/projects/3.CIIR/2.benchmark/run/inferCNV/merged_data/gencode.v32.annotation.gtf.RData" )

GTFfile="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/inferCNV/merged_data/gencode.v32.annotation.gtf.RData"
load(GTFfile)

grch38v32_genes = GTF %>% 
  dplyr::filter( type == "gene") %>%
  dplyr::select(seqnames, start, end, gene_type, gene_name) 
print("How many genes in GTF: ")
print(dim(grch38v32_genes)[1])

grch38v32_genes = grch38v32_genes[!duplicated(grch38v32_genes$gene_name),]
print("How many genes left in GTF: ")
print(dim(grch38v32_genes)[1])

rm(GTF)

print("Load anno_genes")
## make anno for gene
anno_genes = data.frame(
  chr = grch38v32_genes$seqnames,
  start = grch38v32_genes$start,
  end = grch38v32_genes$end,
  row.names = grch38v32_genes$gene_name)

colnames(anno_genes) = NULL

## check genes covered by gtf
## base::setdiff(rownames(counts),grch38v32_genes$gene_name  )

#### spots annotation ####

# annotations_file = data.frame(bcsa1@active.ident)

ref.table = read.table(file=paste( path, "BCSA1_ref_spot.txt", sep = "/"))
tumor.keep = read.table(file=paste( path, "BCSA1_filtered_tumor_70_percent_spots_rep1_for_inferCNV.txt", sep = "/"))

tumor.dotsid = intersect(names(bcsa1@active.ident),tumor.keep$V1)

ref.dots = data.frame(group= "Normal", rownames = ref.table$V1)
tumor.dots = data.frame( group= "Tumor", rownames = tumor.dotsid)

annotations_file = rbind(ref.dots, tumor.dots)
rownames(annotations_file) = annotations_file$rownames

annotations_file = annotations_file %>% dplyr::filter( rownames(.) %in% colnames(bcsa1) ) %>% dplyr::select( group)
colnames(annotations_file) = NULL

#### make infercnv object ####


## get counts for each slice and bind by column
print("bind counts")
layers = c("counts.1.1.1", "counts.2.1.2", "counts.1.2.3", "counts.2.2.4")
layersNames = c("BCSA1TumA1","BCSA1TumA2", "BCSA1TumB1", "BCSA1TumB2" )
names(layersNames) = layers

## subsetting
layersNames = layersNames[grepl("1$", layersNames)]

counts = NULL
for (layer in names(layersNames) ){
  
  counts.temp = SeuratObject::LayerData(object =bcsa1[["Spatial"]], layer = layer )
  counts = cbind(counts, counts.temp)
}

dim(counts)

## subsetting
intersect_spots = intersect(colnames(counts),rownames(annotations_file) )
counts = counts[,intersect_spots]
annotations_file = annotations_file[ intersect_spots ,]

## removing HLA
counts = counts[ !grepl("HLA-",rownames(counts) ),]

dim(counts)

BCSA_infCNV = infercnv::CreateInfercnvObject(raw_counts_matrix=counts,
                                             gene_order_file=anno_genes, 
                                             annotations_file=annotations_file, 
                                             ref_group_names=c("Normal"), 
                                             chr_exclude = c("chrM"))



BSCAS_infCNV_out = infercnv::run(BCSA_infCNV,cutoff=0.1, 
                                 out_dir= paste(respath,"/","BCSA1_i3", sep = "" ), 
                                 cluster_by_groups=TRUE, #We want to run the analysis by spot, so need this parameter = TRUE
                                 denoise=TRUE,
                                 HMM=TRUE, #We need the HMM data for the visualization
                                 analysis_mode = "subcluster", #We want to run the analysis by spot (inferCNV was designed for scRNAseq)
                                 HMM_report_by = "cell",
                                 HMM_type = "i3",
                                 num_threads = 8,
                                 output_format = "pdf")
  







