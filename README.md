## BreastCancer CTA

Computational tissue annotation (CTA) can provide high-resolution annotation on the H&E images paired with spatial transcriptomics data in breast cancer samples. In this vignette we will use `CTA_align()` function to align CTA output to the corresponding Visium spots. 

### 1. Load the dependencies and souce the function script

We start by loading the packages needed for the analyses. Please install them if you haven't.

```
library(Seurat)
library(tidyverse)
library(data.table)

# Source the function script
source("~/CTA_align.R")

```
### 2. Load the dataset and output from CTA

The example datasets can be found in the example folder. We demonstrate the tutorial using one of the breast cancer samples from patient 2 (BCSA2) in our paper. CTA_output.txt can be obtained after training the classifier using QuPath (**Measure** --> **Show detection measurements** --> **Save**).

```
# Load the Seurat object for the Visium data
load("~/example/ST_data.RData")

# Read the CTA output file
CTA_out <- read.delim("~/example/CTA_output.txt")

```

### 3. Align CTA output to the corresponding Visium spots

In the table below, available parameters in `CTA_align()` are shown:

| Parameters       | Description                                                                                                | 
|------------------|------------------------------------------------------------------------------------------------------------|
| ST_data          | Seurat spatial transcriptomics data object from Visium experiment                                          | 
| image_name       | Name of the image in the ST_data@images slot                                                               | 
| scaling_factor   | Scaling factor used if the raw microscope output images is downsized and used as input for SpaceRanger     | 
| pixel_size       | Pixel size of the tiff images in μm/pixel                                                                  | 
| CTA_output       | A dataframe object containing the detailed information such as classification and coordinates from QuPath 'show detected measurements' function; To optimize reading speed, please use the fread() function from the data.table library   | 
| raw_image_width  | Image width in pixel of the raw microscope output images used for CTA annotation                           | 
| raw_image_height | Image height in pixel of the raw microscope output images used for CTA annotation                          | 
| tiff_image       | Logical value, if tiff image is used for CTA annotation, input TRUE; if jpeg image is used, input FALSE    | 

In our case, we use tiff image to perform the CTA annotation in QuPath, thus, we need to set `tiff_image = TRUE`. We downsize the original image from the microscope to run SpaceRanger, so we use `scaling_factor = 0.3` in the analysis. The pixel size is 0.172 μm/pixel (`pixel_size = 0.1722`). Raw image width is 47616 pixel and raw image height is 48128.

The alignment can be done simply in one function:

```
align_res <- CTA_align(ST_data = bcsa, image_name = "BCSA2TumB1", CTA_output = CTA_out, 
                       scaling_factor = 0.3, pixel_size = 0.172, 
                       raw_image_width = 47616, raw_image_height = 48128, tiff_image = T)
```

The results contains the barcode ID and thus, can be used to, for example, plot the percentage of each cell type on the existing UMAP and calculate the percentage of cell type in each GEX-based clusters. 

