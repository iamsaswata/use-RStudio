library(terra)
library(stars)
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)

#################################################################
## TERRA - https://rspatial.org/terra 
#################################################################

# Classes
# SpatRaster and SpatVector


# Intro

x <- rast()
x

# Create a raster
x <- rast(ncol=36, nrow=18, xmin=-1000, xmax=1000, ymin=-100, ymax=900)
res(x)

# Change the spatial resolution of an existing object
res(x) <- 100
res(x)
ncol(x)

# Change number of columns
ncol(x) <- 18
ncol(x)
res(x)

# Set CRS
crs(x) <- "+proj=utm +zone=48 +datum=WGS84"
x


# Create and fill with values

r <- rast(ncol=10, nrow=10, xmin=0, xmax=100, ymin=0, ymax=200)
values(r) <- 1:ncell(r)
plot(r)

values(r)
hasValues(r)


# check dimension and bbox
dim(r)
xmax(r)
xmin(r)




## Load tiff file
filename <- system.file("ex/meuse.tif", package="terra")
r <- rast(filename)

# create three identical SpatRaster objects
r1 <- r2 <- r3 <- rast(nrow=10, ncol=10)
# Assign random cell values
values(r1) <- runif(ncell(r1))
values(r2) <- runif(ncell(r2))
values(r3) <- runif(ncell(r3))
# Combine
s <- c(r1, r2, r3)
s
names(s)

# Subset a layer
subset(s, 'lyr.1')
s %>% subset(3)



#################################################################
## Raster Algebra
#################################################################
r1 <- rast(ncol=10, nrow=10)
values(r1) <- 1:ncell(r1)
s1 <- sqrt(r1)
s1
s2 <- s1 * r1 + 5
values(r1) <- runif(ncell(r1))
r1 <- round(r1)
r1 <- r1 == 1
values(r1)


# replacement function
s1[r1] <- -0.5
s1[!r1] <- 5
s1[s1 == 5] <- 15
values(s1)


# Multiplication with multi-layer objects with different numbers or layers
# The ‘shorter’ objects are ‘recycled’
# For example, if you multiply a 4-layer object (a1, a2, a3, a4) with a 2-layer object (b1, b2):
# the result is a four-layer object (a1b1, a2b2, a3b1, a3b2).

r <- rast(ncol=5, nrow=5)
r
values(r) <- 1
s <- c(r, r+1)
q <- c(r, r+2, r+4, r+6)
x <- r + s + q
x


## Summary functions (min, max, mean, prod, sum, Median, cv, range, any, all) always return a SpatRaster object
a <- mean(r, s, 10)
a
b <- sum(r, s)
b
st <- c(r, s, a, b)
st
sst <- sum(st)
sst

## Use global if instead of a SpatRaster you want a single number summarizing the cell values of each layer.

global(st, 'sum')
global(sst, 'sum')




#################################################################
## High-level methods
#################################################################

# Modifying a SpatRaster object
# func -> crop, trim, merge, aggregate, disagg, resample, shift, extent, warp



# 
r <- rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
values(r) <- 1:ncell(r)
res(r)
ra <- aggregate(r, 2)
res(ra)
r1 <- crop(r, ext(0, 5, 0, 5))
r2 <- crop(r, ext(4, 10, 4, 10))
m <- merge(r1, r2, filename='test.tif', overwrite=TRUE)
plot(m)
plot(r1)
plot(r2)
plot(r)
plot(ra)

# bf, rotate, t

plot(ra)
plot(rotate(ra))
ra
rotate(ra)

# app
r <- rast(ncol=3, nrow=2)
values(r) <- 1:ncell(r)
values(r)

s <- app(r, fun=function(x){ x[x < 4] <- NA; return(x)} )
as.matrix(s)

# lapp
t <- lapp(c(r, s), fun=function(x, y){ x / (2 * sqrt(y)) + 5 } )
as.matrix(t)

# Mask
u <- mask(r, t)
as.matrix(u)

# == 
v <- u==s
as.matrix(v)

# cover
w <- cover(t, r)
as.matrix(w)

x <- classify(w, c(0,2,1,  2,5,2, 4,10,3))
as.matrix(x)


print("1")
#

























#################################################################
## http://www.wvview.org/ossa/ossa/15_Raster_Analysis_terra.html
#################################################################

## Build raster from scratch


basic_raster <- rast(ncols = 3, nrows = 3, xmin = 10, xmax = 12, ymin = 20, ymax = 22, crs="+init=EPSG:4326")
vals <- seq(1, 9) #sample(c(0:255), 100, replace=TRUE)

values(basic_raster) <- vals

tm_shape(basic_raster)+
  tm_raster() +  tm_layout(legend.outside = TRUE)

