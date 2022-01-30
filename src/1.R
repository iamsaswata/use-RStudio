library(terra)
library(stars)
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)

#################################################################
## http://www.wvview.org/ossa/ossa/15_Raster_Analysis_terra.html
#################################################################

## Build raster from scratch


basic_raster <- rast(ncols = 3, nrows = 3, xmin = 10, xmax = 12, ymin = 20, ymax = 22, crs="+init=EPSG:4326")
vals <- seq(1, 9) #sample(c(0:255), 100, replace=TRUE)

values(basic_raster) <- vals

tm_shape(basic_raster)+
  tm_raster() +  tm_layout(legend.outside = TRUE)

