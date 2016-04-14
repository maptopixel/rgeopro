# abstract = Reclassifies raster based on a specificed number of quantile class breaks.
# Inp: Real-valued raster
#      Desired number of quantiles
# Out: Raster classified by quantiles
# Author: Julian Rosser
###############################################################################

# wps.des: id = geo.reclass-quantile, title = Takes a raster and reclassifies based on quantiles,
# abstract = Takes a real valued raster and reclassifies based on a specificed number of quantile class breaks;
# wps.in: inputRasterModel, type = tiff;
# wps.in: nClassBreaks, double, value=2;

###############################################################################

library(raster); library(classInt);

#Read in inputs from the exec request, this must come before any wd changes!
inputRaster = readGDAL(inputRasterModel)

initialWd = getwd()

# wps.off;

#hardcoded test variables
inputRasterName = paste0("C:/Users/ezzjfr/post_doc/data/oxford_lidar/DTM","/", "ox_Lidar_AND_osDTM5_5m_esri_slope.tif")
inputRasterName = paste0("C:/Users/ezzjfr/post_doc/data/oxford_flickr/cumul_views/", "cumulativeViewshed_dtm5_NEW_500m_1st_Jan_upto_11th_Jan.tif")
#inputRasterName = paste0("C:/temp3/","flickInter.tif")

inputRaster = readGDAL(inputRasterName)
nClassBreaks = 6
# finished with hardcoding test variables, turn wps back on

# wps.on;

#Calculate quantiles
r = raster(inputRaster)
rValues = as.numeric(r@data@values)
#rValues = as.integer(as.numeric(r@data@values))

#remove all the 0 valued
removedRValues = rValues[rValues > 0] #bit of a hack needed for cumulative viewshed
classBreaks = classIntervals(removedRValues, n=nClassBreaks, style="quantile")
classBreaks$brks

#binQ = bins.quantiles(rValues, nClassBreaks,12)
#binQ$binct
#table(cut2(rValues, g=nClassBreaks,levels.mean = TRUE))    # quantile groups


classBreaksVector = classBreaks$brks
classBreaksVector[1] = -Inf
classBreaksVector[length(classBreaksVector)] = Inf

classBreaksVectorCol1 = classBreaksVector[1:length(classBreaksVector)-1]
classBreaksVectorCol2 = classBreaksVector[2:length(classBreaksVector)]
classBreaksVectorCol3 = 1:(length(classBreaksVector)-1)

rclmat  = cbind(classBreaksVectorCol1, classBreaksVectorCol2, classBreaksVectorCol3)
rclmat = data.matrix(rclmat)
rclmat
rc <- reclassify(r, rclmat)
rc
#table(rc@data@values)


#export data to file
outFile = "reclassRasterSlope.tif";
outputRasterModel = writeGDAL(as(rc, 'SpatialGridDataFrame') , outFile, drivername="GTiff", type="Int32")

# wps.out: id = outputRasterModel, type = tiff, title = Reclassifier raster;

