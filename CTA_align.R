#' Align single-cell level annotation from computational tissue annotation (CTA) to Visium spots
#' 
#' @import Seurat
#' @import tidyverse
#' @import data.table
#' 
#' @param ST_data Seurat spatial transcriptomics data object from Visium experiment
#' @param image_name Name of the image in the ST_data@images slot
#' @param scaling_factor Scaling factor used if the raw microscope output images is downsized and used as input for SpaceRanger
#' @param pixel_size Pixel size of the tiff images in μm/pixel
#' @param CTA_output A dataframe object containing the detailed information such as classification and coordinates from QuPath 'show detected measurements' function; To optimize reading speed, please use the fread() function from the data.table library
#' @param raw_image_width Image width in pixel of the raw microscope output images used for CTA annotation
#' @param raw_image_height Image height in pixel of the raw microscope output images used for CTA annotation
#' @param tiff_image Logical value, if tiff image is used for CTA annotation, input TRUE; if jpeg image is used, input FALSE
#' 
#' @return Data frame containing cell counts and percentage for each Visium spot
#' 
#' @author Tianyi Li
#' @references reference
#' 
#' @examples 
#' 

CTA_align <- function(ST_data, image_name, scaling_factor, pixel_size, CTA_output, 
                      raw_image_width, raw_image_height, tiff_image) {
  
  # Get image coordinates from the Seurat object
  df <- get(image_name, ST_data@images)
  temp1.coor <- df@coordinates
  temp1.coor$coor <- paste0(temp1.coor$imagerow, "x", temp1.coor$imagecol)
  temp1.coor$cellid <- rownames(temp1.coor)
  
  # Get the width and height for images used as input for SpaceRanger
  #width_sr <- temp1.coor$imagecol
  #height_sr <- temp1.coor$imagerow
  
  # Get the width and height for raw microscope images used for CTA annotation
  width <- raw_image_width
  height <- raw_image_height
  
  # Calculate the full resolution coordinates, if the microscope images are used for SpaceRanger, the scaling factor is 1
  scaling_factor <- scaling_factor
  temp1.coor$fullcol <- temp1.coor$imagecol / scaling_factor
  temp1.coor$fullrow <- temp1.coor$imagerow / scaling_factor
  
  # Calculate the radius of the Visium spot in pixel
  r <- df@scale.factors$spot
  radius <- (r / 2) / scaling_factor
  
  # Reverse the y coordinates as Visium output normally has y coordinated reverted
  temp1.coor$fullrow_reverse <- height - temp1.coor$fullrow
  
  # To align the coordinates between CTA output and Visium data, the coordinates need to be in pixel
  # Coordinates of tiff images (in um) need to be converted to pixel using the pixel size
  if (tiff_image == TRUE) {
    CTA_output$x <- CTA_output$Centroid.X.µm / pixel_size
    CTA_output$y <- CTA_output$Centroid.Y.µm / pixel_size
  } else {
    CTA_output$x <- CTA_output$Centroid.X.px
    CTA_output$y <- CTA_output$Centroid.Y.px
  }
  
  # Reverse the y axis to oriantate the image
  CTA_output$y_reverse <- height - CTA_output$y
    
  # Calculate the block coordinate for row and column, row is y axis, col is x axis
  temp1.coor$ymin <- temp1.coor$fullrow_reverse - radius
  temp1.coor$ymax <- temp1.coor$fullrow_reverse + radius
  temp1.coor$xmin <- temp1.coor$fullcol - radius
  temp1.coor$xmax <- temp1.coor$fullcol + radius
  
  # If the cell falls in the square of the Visium capture area, then assign the cell id to the annotated cells
  CTA_output$cellid <- NA
  for (j in 1:nrow(CTA_output)) {
    x <- CTA_output$x[j]
    y <- CTA_output$y_reverse[j]
    temp4 <- temp1.coor[which(x > temp1.coor$xmin & x < temp1.coor$xmax & y > temp1.coor$ymin & y < temp1.coor$ymax), ]
    if (nrow(temp4) == 0) {
      CTA_output$cellid[j] <- "Not_captured"
    } else {
      CTA_output$cellid[j] <- temp4$cellid
    }
  }
  
  # Remove the cells not captured
  d4 <- CTA_output %>% filter(cellid != "Not_captured")
  
  # Count number of tumor, immune, and stroma cells from each Visium spot
  d5 <- d4 %>% group_by(cellid) %>% 
    dplyr::select(cellid, Classification, x, y_reverse) %>% 
    dplyr::count(Classification) %>% 
    dplyr::filter(Classification != "") %>% 
    pivot_wider(names_from = Classification, values_from = n)
  d5[is.na(d5)] <- 0
  
  # Change colnames from Immune Cells to Immune
  colnames(d5)[which(colnames(d5) == "Immune cells")] <- "Immune"
  
  # Calculate the total cell count and percentage
  d5$total <- rowSums(d5[,-1])
  d5$tumor_per <- d5$Tumor / d5$total
  d5$immune_per <- d5$Immune / d5$total
  d5$stroma_per <- d5$Stroma / d5$total
  
  return (d5)
}
