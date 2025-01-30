#### loading PKGs ####

library(spacexr)
library(Matrix)
library(doParallel)

library(stringr)

library(Matrix)
library(optparse)

#### paths ####

inpath="/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets/PAM50Bulk/results"
refpath="/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets/Wu_etal_2021_BRCA_scRNASeq/stereoscope"
infopath="/home/qiao.yang/projects/3.CIIR/2.benchmark/datasets"
respath="/home/qiao.yang/projects/3.CIIR/2.benchmark/run/RCTD/Results_RCTD"


#### make options ####


option_list = list(
  make_option(c("--patientid"), default=FALSE,type='character',
              help="Which patinetid")
)

opt = parse_args(OptionParser(option_list=option_list))

patientid = opt$patientid


#### Loading info table ####
info = read.table(file = paste0(infopath, "/" ,"BCSA_info.txt"), sep = "\t", header = T )


#### run RCTD #### 

## loading reference
if (info$type[ which(info$id == patientid) ] == "TNBC") {
  load(file = paste( refpath,"tnbc", "ref_RCTD_tnbc.RData", sep = "/"))
} else if (info$type[ which(info$id == patientid) ] == "HER2+") {
  load(file = paste( refpath,"her2", "ref_RCTD_her2.RData", sep = "/"))
}

## loading data
load(paste0(inpath, "/", patientid, ".RData"))

## filter bsca and get matrix
bcsa = subset(bcsa, subset = nFeature_Spatial > 500 & percent_mito < 25 & percent_hb < 0.2)

mtx_input = bcsa[["Spatial"]]$counts

## coordinates
coordsmtx = get(patientid, bcsa@images)
coordsmtx = coordsmtx@coordinates[,c("row","col")]


##UMI
UMImtx = colSums(mtx_input)

## object
mtx = SpatialRNA(coordsmtx, mtx_input, UMImtx)


## RCTD
myRCTD = create.RCTD(mtx, reference, max_cores = 4,CELL_MIN_INSTANCE = 2) 
myRCTD = run.RCTD(myRCTD, doublet_mode = 'full')

barcodes = colnames(myRCTD@spatialRNA@counts)
weights = myRCTD@results$weights
norm_weights = normalize_weights(weights)
norm_weights = as.matrix(norm_weights)

write.table(norm_weights, file =  paste0( respath, "/",patientid,"/", patientid ,"_norm.txt"),
            quote = F, sep = "\t", row.names = T)

save(myRCTD,file = paste0( respath, "/",patientid,"/", patientid ,"_full.RData") )


