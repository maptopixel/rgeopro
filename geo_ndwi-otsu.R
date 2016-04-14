# abstract = Takes a real valued raster bands and computes the NDWI and then thresholds the index map;
# Inp: Real-valued raster stack
#      Threshold (optional)
# Out: Raster classified by quantiles
# Author: Julian Rosser
###############################################################################

# wps.des: id = geo.ndwi-otsu,
# title = takes raster bands for computing NDWI,
# abstract = Takes a real valued raster bands and computes the NDWI and then thresholds the index map;
# wps.in: inputRasterModel, type = tiff;

###############################################################################

library(raster); library(classInt); library(rgdal)


#Read in inputs from the exec request, this must come before any wd changes!
inputRaster = readGDAL(inputRasterModel)

initialWd = getwd()

# wps.off;

#hardcoded test variables
inputRasterName = paste0("C:/Users/ezzjfr/post_doc/data/spatial_accuracy","/", "lsat_DOSReflectance_b1-7_spatial_accuracy.tif")

inputRasterStack = stack(inputRasterName)
# finished with hardcoding test variables, turn wps back on

# wps.on;

#rf <- writeRaster(inputRasterStack, filename="test.tif", format="GTiff", overwrite=TRUE)

#Calculate NDWI
b3 <- raster(inputRasterStack, layer=3)
b4 <- raster(inputRasterStack, layer=4)
b5 <- raster(inputRasterStack, layer=5)
b7 <- raster(inputRasterStack, layer=7)


#ndviRaster = (b5-b4) / (b5 + b4) # NDVI
mndwiRaster = (b3-b7) / (b3 + b7) # Modified ndwi (better for openwater)


#Could use something like -0.04 or 0.04 depending on numerator OR use otsu
#otsuThreshold = calculateOtsu(mndwiRaster@data@values) #function broken
#otsuThreshold


#ndwiRasterOtsu = ndwiRaster < otsuThreshold


otsuThreshold = 0.2667 #matlab in-built Otsu calculated this threshold using the full (i.e. 6174000 pixel area) MNDWI image
#otsuThreshold = 0.6278945 #CBimage R package calculated this


# set values below 100 to NA.
fun <- function(x) { x[x>otsuThreshold] <- 1; x[x<=otsuThreshold] <- 0; return(x) }
rc2 <- calc(mndwiRaster, fun)



#export the data
outFile = "mndwiRasterDOSReflectance.tif";
#export the actual mndwi layer
outputRasterModel = writeGDAL(as(mndwiRaster, 'SpatialGridDataFrame') , outFile, drivername="GTiff")


outFile = "mndwiRasterDOSReflectanceNDVI_otsu_matlab.tif";
#export the thresholded value 
outputRasterModel = writeGDAL(as(rc2, 'SpatialGridDataFrame') , outFile, drivername="GTiff") # this should probably be an integer raster

# wps.out: id = outputRasterModel, type = tiff, title = Reclassifier raster;

