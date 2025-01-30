library(Giotto)
library(tidyverse)
library(Seurat)


respath="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/spDWLS/res_spDWLS/"


# Load meta data
info <- openxlsx::read.xlsx("/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets/ST-AR_Visium_samples_with file names.xlsx") %>% na.omit()
info$id <- paste0(info$Patientid, info$Tumor.area)
info <- info[-which(info$id == "BCSA2TumA1"), ]
info <- info[-which(info$id == "BCSA2TumA2"), ]

# Create scRNAseq reference for TNBC samples
matrix <- read_tsv("/Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/tnbc/latest/BCSA1TumA1_tnbc_mtx.tsv") %>% as.data.frame()
meta <- read_tsv("/Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/tnbc/latest/BCSA1TumA1_tnbc_meta.tsv")

rownames(matrix) <- matrix[,1]
matrix <- matrix[, -1]
tnbc_matrix <- t(matrix) %>% as.data.frame()

# Read highly variable genes and cell type markers for deconvolution
feature <- read.delim("/Users/tili/Desktop/CIIR/data/ST_data/stereoscope/features.txt", header = F)
#tnbc_matrix <- tnbc_matrix[which(rownames(tnbc_matrix) %in% feature$V1),]

giotto_SC <- createGiottoObject(
  expression = tnbc_matrix,
  instructions = NULL
)

giotto_SC <- addCellMetadata(giotto_SC,
                             new_metadata = meta)

giotto_SC <- normalizeGiotto(giotto_SC)

tnbc_DWLS_matrix <- makeSignMatrixDWLSfromMatrix(matrix = getExpression(giotto_SC, values = "normalized", output = "matrix"),
                                                 cell_type = pDataDT(giotto_SC)$bio_celltype,
                                                 sign_gene = feature$V1)

# Create scRNAseq reference for HER2 samples
matrix <- read_tsv("/Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/her2/latest/BCSA2TumB1_her2_mtx.tsv") %>% as.data.frame()
meta <- read_tsv("/Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/her2/latest/BCSA2TumB1_her2_meta.tsv")

rownames(matrix) <- matrix[,1]
matrix <- matrix[, -1]
her2_matrix <- t(matrix) %>% as.data.frame()

# Keep only the genes used for deconvolution
feature <- read.delim("/Users/tili/Desktop/CIIR/data/ST_data/stereoscope/features.txt", header = F)

giotto_SC <- createGiottoObject(
  expression = her2_matrix,
  instructions = NULL
)

giotto_SC <- addCellMetadata(giotto_SC,
                             new_metadata = meta)

giotto_SC <- normalizeGiotto(giotto_SC)

her2_DWLS_matrix <- makeSignMatrixDWLSfromMatrix(matrix = getExpression(giotto_SC, values = "normalized", output = "matrix"),
                                                 cell_type = pDataDT(giotto_SC)$bio_celltype,
                                                 sign_gene = feature$V1)


#################################### Perform deconvolution

## former part 1:23
for (i in 1:23) {
  # Load seurat ST data
  # The Processed data have <25% mitochondrial reads, <20% hb-reads, at least 500 detected genes, remove mitochondria genes
  patientid <- info$id[i]
  print(patientid)
  load(paste0("/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets/PAM50Bulk/results/", patientid, ".RData"))
  
  # Filter the low quality spots and find highly variable genes
  bcsa <- subset(bcsa, subset = nFeature_Spatial > 500 & percent_mito < 25 & percent_hb < 0.2)
  
  # set working directory
  results_folder = "/home/qiao.yang/projects/3.CIIR/2.benchmark/run/spDWLS/res_spDWLS"
  
  # Optional: Specify a path to a Python executable within a conda or miniconda
  # environment. If set to NULL (default), the Python executable within the previously
  # installed Giotto environment will be used.
  my_python_path = NULL # alternatively, "/local/python/path/python" if desired.
  
  # Create Giotto Instructions
  instrs = createGiottoInstructions(save_dir = file.path(results_folder, patientid, "giotto.output"),
                                    save_plot = TRUE,
                                    show_plot = TRUE,
                                    python_path = my_python_path)
  
  ## provide path to visium folder
  data_path = "/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets/visium_spatial_out/"
  
  ## directly from visium folder
  gbcsa <- createGiottoVisiumObject(visium_dir = file.path(data_path, info$Visium.ID[i],"/outs"),
                                    expr_data = 'raw',
                                    png_name = 'tissue_lowres_image.png',
                                    gene_column_index = 2,
                                    instructions = instrs)
  
  # Keep cells that pass the filtering
  cellid <- colnames(bcsa)
  gbcsa <- subsetGiotto(gbcsa, cell_ids = cellid)
  
  ## normalize
  gbcsa <- normalizeGiotto(gobject = gbcsa, scalefactor = 10000, verbose = T)
  
  ## add gene & cell statistics
  gbcsa <- addStatistics(gobject = gbcsa)
  
  ## highly variable features / genes (HVF)
  gbcsa <- calculateHVF(gobject = gbcsa, save_plot = F)
  
  cluster <- bcsa@meta.data
  gbcsa <- addCellMetadata(gbcsa,
                           new_metadata = cluster)
  
  if (info$type[i] == "TNBC") {

    gbcsa <- runDWLSDeconv(gobject = gbcsa, sign_matrix = tnbc_DWLS_matrix, cluster_column = "spaceranger_clusters", n_cell = 50)
  } else {
    
    gbcsa <- runDWLSDeconv(gobject = gbcsa, sign_matrix = her2_DWLS_matrix, cluster_column = "spaceranger_clusters", n_cell = 20)
  }
  
  dwls <- gbcsa[["spatial_enrichment"]][[1]]@enrichDT
  write.csv(dwls, file = paste0(results_folder,"/", patientid, "_SpatialDWLS.csv"))
  
}
