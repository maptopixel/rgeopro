# Uncertainty/sensitivity analysis for Spatial Accuracy 2016
# Leibovici DG & Rosser JF  (2016) Multiway sensitivity analysis of the fusion of earth observation, topography and social media data for rapid flood mapping. Spatial Accuracy, 12th International Symposium, 5-8 July 2016, Montpellier, France
#  of the workflow for rapid inundation extent estimation
#   from  Rosser JF, Leibovici DG, Jackson MJ (2016) Rapid flood inundation mapping using social media, remote sensing and topographic data.  Natural Hazards (submitted)
#   
# Dr Didier Leibovici, Dr Julian Rosser and Pr Mike Jackson University of Nottingham
#
#  estimate/use fixed uncertainty  
#     then perform the sampling 
#     then the Woe to estimate the posterior prob
#     then compute the summaries of the n simulations per pixel as
#     mean, variance and quartiles (Q1 Median Q3)
#   to be used a 5 way array for multiway sensitivity analysis
#######################################################################
# nohup nice  -n -15 R CMD BATCH UncertaintySensitivity.R &
#########

library(rgeos); library(maptools); library(rgdal); 
library(spatstat); library(rgrass7); library(raster); library(sp);

#setwd("/Users/lgzdl/Documents/Dids/DidsE/R_GIS/R_SpatialAccuracy2016/")

setwd("./DataSpatialAccuracy")
# uncertainty dim
#
# 1st dim 	Social flickr 
       # point resampled need the whole raster for viewshed
# 2nd dim 	Topo slope and elevation  slope and elevation from initial DEM then attributes resampled
# 3rd dim   EO  MNDWI attribute resampled
  #########
# 4th dim will be measures
# 5th dim will be map/raster 


# defined classes cuts from the initial fitting of the WoE
load("rclmat.flickR.RData")
rclmat.flickR=rclmat
load("rclmat.mndwi.RData")
rclmat.mndwi=rclmat
load("rclmat.slope.RData")
rclmat.slope=rclmat
load("rclmat.elevation.RData")
rclmat.elev=rclmat
rm(rclmat)
# read data and set/compute uncertainties levels 1vr 2vr 3vr 4vr for the 3 dimensions
 
 # flickR (note we use the the same DTM here)
 
#Global variables for type of processing
outputRasterIsViewshedCount = 1; #0 = simple logical union of viewsheds. 1 = cumulative viewshed (viewshed count)
maxDistance = 500; #maximum viewing distance in metres for each geo-point.
writeViewshedsToDisk = 0; #for debugging after each viewshed calcuation by GRASS

 ##
	flickr="Flickr_for_spatial_accuracy"
    ViewPoints = readOGR(dsn=paste0(flickr, ".shp"), layer= flickr)
	
	baseDemFilename="Lidar_DTM_downsampled_to_5m_spatial_accuracy.tif" 			
    baseElevationModel = readGDAL(baseDemFilename)
    demProjection = proj4string(baseElevationModel)
#/Applications/GRASS-6.4.app/Contents/MacOS/bin/
#########################
loc <- initGRASS("/usr/lib/grass70",home=getwd(), gisDbase="GRASS_TEMP", override=TRUE )

#grasspath="/Applications/GRASS-6.4.app/Contents/MacOS/bin/"

execGRASS("r.in.gdal", flags="o", parameters=list(input=baseDemFilename, output="DEM"))

execGRASS("g.region", parameters=list(raster="DEM"))

 
 # grid
grd <- gmeta2grd()
ncells <- grd@cells.dim[1]*grd@cells.dim[2] #759980

bottomLeftX = grd@cellcentre.offset[1]
bottomLeftY = grd@cellcentre.offset[2]
###########################################

#demBbox = bbox(baseElevationModel)
BasedemExtentPoly = as(extent(baseElevationModel), "SpatialPolygons")
###
flickR.u=c(2.5,5,10,15) # as 68% halo (i.e 1sd) 
###
 
CumViewShed.u<-function(VP=ViewPoints,u=flickR.u[1],vect=TRUE,cl=rclmat.flickR, BasedemExtentPolyc= BasedemExtentPoly,ncellsc=ncells,macDistance=maxDistance,grdc=grd, demProjectionc= demProjection){
	#
	mapCalc <- function(los, productRaster){
  #function to add up the viewshed layers
  #these conditions seem to be needed for lidar data
  if (any(is.na(los[[1]])) == TRUE) {
    losValues <- ifelse(is.na(los@data[[1]]), 0, 1) #if na values, set to 0, else set to 1
  } else {
    losValues = los@data[[1]]
  }
  sum <- losValues + productRaster
  return(sum)
}
    #
    rcl<-function(v,cl){
    	#cl with three columns 1=[a 2=,b] or b finite or inf[  3=class nb
    	# so nb rows in cl are the nb of classes
    	v=as.vector(v)
    	for(c in 1:dim(cl)[1]) v[ v>cl[c,1] & v<=cl[c,2] ]=cl[c,3]
    return(v)	
    }
#sample

	AA = coordinates(VP) +matrix(rnorm(prod(dim(VP@coords)),0,u),nrow=dim(VP@coords)[1], ncol=2)
	VP=SpatialPoints(AA, proj4string=CRS(proj4string(VP)), bbox = NULL)
	
	#crop the points to the dem extent
demExtentPoly = BasedemExtentPolyc#as(extent(baseElevationModel), "SpatialPolygons") 
	 proj4string(demExtentPoly) = proj4string(VP)
    VP = crop(VP,demExtentPoly)

#ViewPoints prep
	aFrame = data.frame(VP) # create data frame
#coords = aFrame[2:3] # take the coordinate cols
	coordsX = aFrame[["coords.x1"]]
	coordsY = aFrame[["coords.x2"]]
	coords = cbind(coordsX,coordsY)
#  ini	 
sumV <- rep(0, ncellsc)
uniqueOwners = unique(aFrame$ownerNum)

#First viewshed
execGRASS("r.viewshed", parameters = list(input = "DEM", output = "cumulativeViewshed", max_distance= macDistance, coordinates = as.integer(coords[1,])), flags = c("overwrite" , "b","quiet"))

los <- readRAST("cumulativeViewshed")
#add this firest layer into the cumulative raster
sumV = mapCalc(los,sumV)

# Now loop over any remaining points, creating viewsheds and adding into cumulative viewshed with logical union mapcalc

for (i in 2:nrow(coords)) {
#for (i in 2:4) {
  losLayerName = "cumulativeViewshed"#toString(i)
  #initial cleaning
  #execGRASS("g.remove", parameters = list(type = "raster",name=losLayerName), flags="f")
  #rm(losInR)

  #Compute this shed instance
  execGRASS("r.viewshed", parameters = list(input = "DEM", output = losLayerName, max_distance=maxDistance, coordinates = as.integer(coords[i,])), flags = c("overwrite" , "b","quiet"))

  los <- readRAST(losLayerName)

    #add it to the cumulative viewshed with map calc
  sumV = mapCalc(los,sumV)
  
}#end of for
#reclass 
sumV=rcl(sumV,cl)

# map to sp	
	if(!vect){
	sgdf <- SpatialGridDataFrame(grdc, data = data.frame(sum=sumV))
	proj4string(sgdf) = demProjectionc
		return(sgdf)
	}
	else{
		 return(sumV)
		}
} #CumViewShed.u
 

 # slope and elevation and MNDWI
   slope.ini = as.vector(as.matrix(readGDAL("Lidar_DTM_slope_spatial_accuracy.tif")@data))
   slope.ini[is.na(slope.ini)]=0 # NA
  elev.ini= as.vector(as.matrix(baseElevationModel@data))
         lsat=readGDAL("lsat_DOSReflectance_b1-7_spatial_accuracy.tif")
 #ndviRaster = (b5-b4) / (b5 + b4) # NDVI
 #mndwiRaster = (b3-b7) / (b3 + b7) # Modified ndwi (better for openwater)
  mndwi.ini=(lsat@data$band3 -lsat@data$band7)/(lsat@data$band3 +lsat@data$band7)
 ###uncertainty global ...sd(slope.ini)=2.98 using the observed sd as base
 slope.u=c(1.5,3,6,9) # 1/2sd, 1sd, 2sd,3sd
 elev.u=c(8,16,32,48)
 mndwi.u=c(0.11,0.22,0.44,0.66)
 ###
 rnclass.vect<-function(vect= slope.ini,u=slope.u[1],cl=rclmat.slope){
 	  rcl<-function(v,cl){
    	#cl with three columns 1=[a 2=,b] or b finite or inf[  3=class nb
    	# so nb rows in cl are the nb of classes
    	for(c in 1:dim(cl)[1]) v[ (v>cl[c,1] & v<=cl[c,2]) ]=cl[c,3]
    return(v)	
    }#rcl
     vect=as.vector(vect)
 	 vect=vect+rnorm(length(vect),0,u)
 	 return(rcl(vect,cl))
 }#end of rnclass
####################### could have been binary for slope and elev
WWflickr=c(-0.23,-0.15,-0.15,1.72,1.6)
WWmndwi=c(-1.57,2.11)
WWslope=c(0.37,0.37,-1.41,-1.41,-1.41,-1.41)
WWelev=c(1,1,-4.15,-4.15,-4.15,-4.15)
WW=list(WWflickr, WWmndwi, WWslope, WWelev)
Watprior=  2e-05  #   nb of traning sample or 0.09 579843/(579843+5306975) mndwi prior? 
#######################100/(579843+5306975) 1.7e-05
WoEp<-function(flickR, mndwi,slope,elev,wei=WW,prior=Watprior){
	# WoE has estimated the weights and here is the calculation of the posterior
	#see git WoeR
	#WW vector of wieghts for the classes but once transformed into binary (could transform it into binary but ...) so it is either the W+ or W-
	 # logit(Water(px)/f,s,e,m)=logit(priorWater)+sum_i(cl(i(px))
	 attW<-function(v,w){
	 	#v as 1:cl where cl is lenght of w
	 	for(i in 1:length(w))v[v==i]=w[i]
	 return(v)
	 }#attW  exp(logit)
	 el=exp(log( prior/(1-prior) )+attW(flickR,wei[[1]])+attW(mndwi,wei[[2]])+attW(slope,wei[[3]])+attW(elev,wei[[4]]))
     posterior=el/(1+el)
return(posterior)
}
##################
# simul (first version without parallelisation)
#nDsimul=100 # per dimension making 1 000 000 for 100 (same nb in each dimension)
#or
nDsimul=100 # evaluations per combination of uncertainties
u1=4;u2=4;u3=4
Sensi=array(0,dim=c(ncells,u1,u2,u3,6),dimnames=list(NULL,paste("soc",1:u1,sep=""),paste("top",1:u2,sep=""),paste("eo",1:u3,sep=""),c("Min","Q1","Med","Mean","Q3","Max")))


###########################


###snow cluster prep
 library(snow)
  #CLusters=rep("localhost",16) # maybe something different ???  type of clusters mpi etc...
   cl= makeSOCKcluster(2)
     clusterEvalQ(cl, library(rgeos));  clusterEvalQ(cl, library(maptools));  clusterEvalQ(cl, library(rgdal));  clusterEvalQ(cl, library(spatstat)); clusterEvalQ(cl, library(rgrass7));
     clusterEvalQ(cl, library(raster)); clusterEvalQ(cl, library(sp));
        
     
 ##########
simul<-function(n,VPs= ViewPoints, fu= flickR.u,socs=soc, rclfR= rclmat.flickR,slo=slope.ini, su= slope.u, tops= top, rcls= rclmat.slope, elevi= elev.ini,elu=elev.u, rcle= rclmat.elev, mndwii= mndwi.ini, mu= mndwi.u, eos= eo, rcli= rclmat.mndwi,weis=WW,priors=Watprior){
		#all objects are looked for the parent frame
		
	# simul flickR then cum viewshed then classes	rclmat
       flickR=CumViewShed.u(VP=VPs,fu[socs],vec=TRUE,cl=rclfR) # vector
	# simul slope and elevation then classes
	   slope=rnclass.vect(vect=slo, u=su[tops],cl=rcls)
	   elev=rnclass.vect(vect=elevi,u=elu[tops],cl= rcle)
	# simul MnDWI then classes
	   mndwi=rnclass.vect(vect=mndwii,u=mu[eos],cl=rcli)
	      # WoE then 
	      #resul[,n]=
	   out=WoEp(flickR, mndwi,slope,elev,wei=weis,prior=priors)   
	      return(out)
	}#simul 
clusterExport(cl,"simul")
    clusterExport(cl, "ViewPoints");clusterExport(cl, "flickR.u");clusterExport(cl, "rclmat.flickR");
     clusterExport(cl, "slope.ini");clusterExport(cl, "slope.u");clusterExport(cl, "rclmat.slope");
    clusterExport(cl, "elev.ini");clusterExport(cl, "elev.u");clusterExport(cl, "rclmat.elev");
     clusterExport(cl, "mndwi.ini");clusterExport(cl, "mndwi.u");clusterExport(cl, "rclmat.mndwi");
clusterExport(cl,"WW"); clusterExport(cl,"Watprior") ; 
clusterExport(cl,"demProjection") ;clusterExport(cl,"grd") ;clusterExport(cl,"maxDistance") ;
clusterExport(cl,"baseElevationModel") ;clusterExport(cl,"BasedemExtentPoly") ;
clusterExport(cl,"ncells") ;clusterExport(cl,"bottomLeftX") ;clusterExport(cl,"bottomLeftY") ;
clusterExport(cl,"writeViewshedsToDisk")
clusterExport(cl,"CumViewShed.u")
clusterExport(cl,"rnclass.vect")
clusterExport(cl,"WoEp")
clusterExport(cl,"baseDemFilename")
#
clusterCall(cl,dir)
#clusterCall(cl, function(){
     	loc <<- initGRASS("/usr/lib/grass70",home=getwd(), gisDbase="GRASS_TEMP", override=TRUE )
     	execGRASS("r.in.gdal", flags="o", parameters=list(input=baseDemFilename, output="DEM"))
     execGRASS("g.region", parameters=list(raster="DEM"))
    # });
#
Toutdeb=date()
cat("debut: ",Toutdeb)

########testing
u1=1;u2=1;u3=1; nDsimul=3
###############
for (soc in 1:u1){ 
for (top in 1:u2){
for (eo  in 1:u3){

	deb=date()
	clusterExport(cl,"soc");clusterExport(cl,"top");clusterExport(cl,"eo");
	#debug(simul)
     resul=parSapply(cl,1:nDsimul,simul)
	cat("simul ",soc,top,eo," :: ") 
	deb ;date()
		
	# summary simul px x 100 0 in x min Q1 Q2 mean Q3 max
   Sensi[,soc,top,eo,]=t(apply(resul,1,summary))
   
}# end of u3
}# end of u2	 
}# end of u1
cat("debut et fin :")
Toutdeb;date()  
  
  SensiGrid=lsat
  #SensiGrid@data=data.frame(resul)
  SensiGrid@data=data.frame(Sensi[,2,2,2,])
 save(SensiGrid, Sensi,file="simulSpAcc2016.RData")
  
# fini pour ça


# mise en grid

