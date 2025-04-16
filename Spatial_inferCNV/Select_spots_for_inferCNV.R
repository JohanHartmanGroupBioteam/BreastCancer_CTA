library(tidyverse)

## Get spots without tumor cells from QuPath for inferCNV
info <- openxlsx::read.xlsx("~/Desktop/CIIR/data/ST-AR_Visium_samples_with file names.xlsx") %>% na.omit()
info$id <- paste0(info$Patientid, info$Tumor.area)
info <- info[-which(info$id == "BCSA2TumA1"), ]
info <- info[-which(info$id == "BCSA2TumA2"), ]

BCSA1_ref <- character()
BCSA2_ref <- character()
BCSA3_ref <- character()
for (i in 1:nrow(info)) {
  patient <- info$Patientid[i]
  id <- info$id[i]
  
  # Read qupath annotation and rename spot ID
  temp1 <- read.csv(paste0("/Users/tili/Desktop/CIIR/results/cell2location/Merged_correlation_results/",
                           id, "_correlation_matrix.csv"))
  temp1$X <- paste0(gsub(patient, "", id), "_", temp1$X)
  ref <- temp1[which(temp1$Tumor == 0),]$X
  
  # Combine the reference spots together
  if (patient == "BCSA1") {
    BCSA1_ref <- c(BCSA1_ref, ref)
  } else if (patient == "BCSA2") {
    BCSA2_ref <- c(BCSA2_ref, ref)
  } else if (patient == "BCSA3") {
    BCSA3_ref <- c(BCSA3_ref, ref)
  }
}

write.table(BCSA1_ref, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA1_ref_spot.txt", 
            row.names = F, col.names = F)
write.table(BCSA2_ref, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA2_ref_spot.txt", 
            row.names = F, col.names = F)
write.table(BCSA3_ref, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA3_ref_spot.txt", 
            row.names = F, col.names = F)

########### Select spots with at least 70% tumor to run inferCNV, remove rep2

# BCSA1
load("~/Desktop/CIIR/results/colocalization/BCSA1_merge.Rdata")
decon <- data.frame()
for (i in names(bcsa1@images)) {
  temp1 <- read.csv(paste0("/Users/tili/Desktop/CIIR/results/cell2location/Merged_correlation_results/", 
                           i, "_correlation_matrix.csv"))
  region <- gsub("BCSA1", "", i)
  temp1$X <- paste0(region, "_", temp1$X)
  decon <- rbind(decon, temp1)
}

bcsa1_cnv_spot <- decon %>% dplyr::filter(tumor_per > 0.7)
bcsa1_cnv_spot$region <- str_sub(bcsa1_cnv_spot$X, end = 5)
bcsa1_cnv_spot <- bcsa1_cnv_spot %>% 
  dplyr::filter(region != "TumA2") %>% 
  dplyr::filter(region != "TumB2")

unique(bcsa1_cnv_spot$region)

write.table(bcsa1_cnv_spot$X, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA1_filtered_tumor_70_percent_spots_rep1_for_inferCNV.txt", col.names = F, row.names = F)

# BCSA2
load("~/Desktop/CIIR/results/colocalization/BCSA2_merge.Rdata")
decon <- data.frame()
for (i in names(bcsa2@images)) {
  temp1 <- read.csv(paste0("/Users/tili/Desktop/CIIR/results/cell2location/Merged_correlation_results/", 
                           i, "_correlation_matrix.csv"))[,c(1:12)]
  region <- gsub("BCSA2", "", i)
  temp1$X <- paste0(region, "_", temp1$X)
  decon <- rbind(decon, temp1)
}

bcsa2_cnv_spot <- decon %>% dplyr::filter(tumor_per > 0.7)
bcsa2_cnv_spot$region <- str_sub(bcsa2_cnv_spot$X, end = 5)
bcsa2_cnv_spot <- bcsa2_cnv_spot %>% 
  dplyr::filter(region != "TumC2") %>% 
  dplyr::filter(region != "TumB2") %>% 
  dplyr::filter(region != "TumD2") %>% 
  dplyr::filter(region != "TumD3") %>% 
  dplyr::filter(region != "TumD4") %>% 
  dplyr::filter(region != "TumE2")

unique(bcsa2_cnv_spot$region)
write.table(bcsa2_cnv_spot$X, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA2_filtered_tumor_70_percent_spots_rep1_for_inferCNV.txt", 
            col.names = F, row.names = F)

# BCSA3
load("~/Desktop/CIIR/results/colocalization/BCSA3_merge.Rdata")
decon <- data.frame()
gro <- names(bcsa3@images)[-4]
for (i in gro) {
  temp1 <- read.csv(paste0("/Users/tili/Desktop/CIIR/results/cell2location/Merged_correlation_results/", 
                           i, "_correlation_matrix.csv"))
  region <- gsub("BCSA3", "", i)
  temp1$X <- paste0(region, "_", temp1$X)
  decon <- rbind(decon, temp1)
}

bcsa3_cnv_spot <- decon %>% dplyr::filter(tumor_per > 0.7)
bcsa3_cnv_spot$region <- str_sub(bcsa3_cnv_spot$X, end = 5)
bcsa3_cnv_spot <- bcsa3_cnv_spot %>% 
  dplyr::filter(region != "TumA2") %>% 
  dplyr::filter(region != "TumB2") %>% 
  dplyr::filter(region != "TumC2") %>% 
  dplyr::filter(region != "TumD2") 
unique(bcsa3_cnv_spot$region)

write.table(bcsa3_cnv_spot$X, file = "/Users/tili/Desktop/CIIR/results/infercnv/BCSA3_filtered_tumor_70_percent_spots_rep1_for_inferCNV.txt", 
            col.names = F, row.names = F)


