#----------------------------------------------------------------------------------------------------
###### GEOG 418 FINAL PROJECT ######

#Load libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(gtable)
library(ggplot2)
library(gridExtra)
library(grid)
library(shinyjs)
library(sp)
library(rgeos)

#Set working directory
dir <- "C:\\Users\\Lucas Aubert\\OneDrive\\University\\Year 4\\Fall 2020\\GEOG 418\\Labs\\Final Project"
setwd(dir)


#----------------------------------------------------------------------------------------------------
###### DATA PREPARATION ######

#Reading in elevation dataset
elev <- readOGR(dsn = './Data', layer = 'ElevSample')
elev <- spTransform(elev, CRS("+init=epsg:26910"))

#Reading in VRI data
VRI <- readOGR(dsn = './Data', layer = 'WatershedVRI')
VRI <- spTransform(VRI, CRS("+init=epsg:26910"))
head(VRI@data)

#clean VRI data and rename columns
vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", "PROJ_AGE_1",
                  "SITE_INDEX", "SPECIES__4", "SPECIES__5",
                  "PROJ_HEI_1", "SPECIES_PC", "SPECIES__6",
                  "VRI_LIVE_S", "BASAL_AREA", "WHOLE_STEM",
                  "CROWN_CL_1")

#Create subset of data
vriClean <- VRI[,vriCleanCols]

#Meta Data (https://www2.gov.bc.ca/assets/gov/farming-natural-resources-and-industry/forestry/stewardship/forest-analysis-inventory/data-management/standards/vegcomp_poly_rank1_data_dictionaryv5_2019.pdf)
# FID = Field ID
# PolyID = VRI Polygon ID
# Stand_Age = Estimated stand age projected to 2020 from estimated establishment date
# Site_Index = A value to estimate site quality. This describes the height that the stand could grow to by age 50 in meters.
# CoDom_Sp = The species code for the co-dominant tree species. Full list of codes: https://www.for.gov.bc.ca/hfp/publications/00026/fs708-14-appendix_d.htm
# Dom_Sp = The species code for the dominant tree species. Full list of codes: https://www.for.gov.bc.ca/hfp/publications/00026/fs708-14-appendix_d.htm
# Stand_HT = The estimated height for the stand
# DomSP_Perc = The estimated percentage of the dominent species
# CDomSP_Perc = The estimated percentage of the co-dominent species
# Stand_Dens = Estimated density of stand (Stems per hectare)
# Stand_BA = Estimated Basal area of the stand (square meters)
# Stand_StemBio = Estimated stand level stem biomass (tonnes per hectare)
# Stand_CrownCl = The percentage of ground area covered by tree crowns

#Rename columns to make them more readable
newNames <- c("FID", "PolyID", "Stand_Age", "Site_Index",
              "CoDom_Sp", "Dom_Sp", "Stand_HT", "DomSP_Perc", 
              "CDomSP_Perc", "Stand_Dens", "Stand_BA", "Stand_StemBio", "Stand_CrownCl")

colnames(vriClean@data) <- newNames
head(vriClean@data)

#Select Site_Index where Site_Index is NOT NA
vriClean@data$Site_Index[vriClean@data$Site_Index == 0] <- NA
vriClean <- vriClean[!is.na(vriClean@data$Site_Index), ]

#Create choropleth map of height
map_SI <- tm_shape(vriClean) +
  tm_polygons(col = "Site_Index",
              title = "Site Index (m)",
              style = "jenks",
              palette = "viridis", n = 5,
              border.alpha = 0.3) +
  tm_shape(elev) + 
  tm_symbols(size = 0.1, col = "red", alpha = 0.6, border.alpha = 0.2) + 
  tm_add_legend(type = "fill", labels = c("Elevation"),
                col = c(adjustcolor("red", alpha.f = 0.6))) +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.35))
map_SI

#tmap colour palettes
tmaptools::palette_explorer()

#----------------------------------------------------------------------------------------------------
###### DESCRIPTIVE STATISTICS ######

##Create a table of descriptive stats

##VRI
#2 most common dominant stand species
mode_Sp = tail(names(sort(table(vriClean$Dom_Sp))), 2)
#1 -> FDC --> Coast Douglas Fir
#2 -> HW --> Western Hemlock

#% of total stands that are dominated by Coast Douglas Fir
num_FDC = length(which(vriClean$Dom_Sp == 'FDC'))
perc_FDC = num_FDC / length(which(vriClean$Dom_Sp != ''))
#-> 80.8% FDC

#% of total stands that are dominated by Western Hemlock
num_HW = length(which(vriClean$Dom_Sp == 'HW'))
perc_HW = num_HW / length(which(vriClean$Dom_Sp != ''))
#-> 7.6% HW

#average stand age
avg.stand.age = round(mean(vriClean@data$Stand_Age), 2)
#average site index
avg.SI = round(mean(vriClean@data$Site_Index), 2)
#range of site index
range.SI = round(range(vriClean@data$Site_Index), 2)
min.SI = range.SI[1]
max.SI = range.SI[2]  

##Elevation
#elev range
range.elev = round(range(elev@data$grid_code), 2)
min.elev = range.elev[1]
max.elev = range.elev[2] 
#avg elev
avg.elev = round(mean(elev@data$grid_code), 2)
#Count elev
num.elev.points = nrow(elev@data)



##Make data tables (Optional)
#Add data to tables
sum_stats_VRI = data.frame(avg.stand.age, min.SI, max.SI, avg.SI)
sum_stats_elev = data.frame(num.elev.points, min.elev, max.elev, avg.elev)


#Make VRItable
VRItable <- tableGrob(sum_stats_VRI) #make a table "Graphical Object" (GrOb) 
VRICaption <- textGrob("Table 1: VRI Summary Statistics (2019)", gp = gpar(fontsize = 10))
padding <- unit(5, "mm")

VRItable <- gtable_add_rows(VRItable, 
                          heights = grobHeight(VRICaption) + padding, 
                          pos = 0)

VRItable <- gtable_add_grob(VRItable,
                            VRICaption, t = 1, l = 2, r = ncol(sum_stats_VRI) + 1)

grid.arrange(VRItable, newpage = TRUE)

#Make elev_table
elev_table <- tableGrob(sum_stats_elev) #make a table "Graphical Object" (GrOb) 
elev_Caption <- textGrob("Table 2: Elevation Summary Statistics (2019)", gp = gpar(fontsize = 10))
padding <- unit(5, "mm")

elev_table <- gtable_add_rows(elev_table, 
                            heights = grobHeight(elev_Caption) + padding, 
                            pos = 0)

elev_table <- gtable_add_grob(elev_table,
                              elev_Caption, t = 1, l = 2, r = ncol(sum_stats_elev) + 1)

grid.arrange(elev_table, newpage = TRUE)

#----------------------------------------------------------------------------------------------------
###### SI SAC ######

##Defining neighbourhood using using Queen's case
vri.nb2 <- poly2nb(vriClean, queen = TRUE)
#Define a network grid
vri.net2 <- nb2lines(vri.nb2, coords=coordinates(vriClean))
#Apply a coordinate system
crs(vri.net2) <- crs(vriClean)

##Create weights matrix
#If it's a neighbour, it will be 1, if not 0
#For example, if there are two neighbours, the weight of 1 will be distributed ovdfer those two - resulting in a weight of 0.5 for each neighbour
#the default weight style is 'w'?
#Queen's case
vri.lw2 <- nb2listw(vri.nb2, zero.policy = TRUE, style = "W")
print.listw(vri.lw2, zero.policy = TRUE)

##Global Moran's i
#Compute Moran's I Test for site index
#Uses the data from weights matrix, polygon values, and zero.policy settings
#Remember, setting zero.policy = True will bypass the calculation of the polygons with no neighbours
#Setting zero.policy = False will result in a crashed program if there are any polygons that have no neighbours in the data (which there are in this case)
mi_SI <- moran.test(vriClean$Site_Index, vri.lw2, zero.policy = TRUE)
mi_SI

#This function will take ~20 minutes to run
#This calculates the expected range of the spatial autocorrelation given connectivity of the polygons (derived fro the weights matrix)
#moran.range <- function(lw) {
#  wmat <- listw2mat(lw)
#  return(range(eigen((wmat + t(wmat))/2)$values))
#}
#moran.range(vri.lw2)


#This extracts the values needed to calculate the z-score
#Global Moran's I
# FOR SOME REASON THIS DOES NOT WORK --> mi_SI <- mi_SI$estimate[[1]]
mi_SI <- 0.4711687
#Expected Moran's I for a random distribution
#eI_SI = mi_SI$estimate[[2]]
eI_SI <- -0.0001790510
#Variance of values in the dataset
#var_SI <- mi_SI$estimate[[3]]
var_SI <- 0.00007331617
#Calculate Z-score
z_SI <- (mi_SI-eI_SI)/var_SI**0.5

#SIGNIFICANT POSITIVE SAC EXISTS!!!! 95% confident

##Make table to show results
#inlcude other SAC results too?

##Local Moran's i

#this will get a local moran's I value for each polygon
lisa.test_SI <- localmoran(vriClean$Site_Index, vri.lw2)
#extracting the results from the lisa test
vriClean$Z.Ii_SI <- lisa.test_SI[,4]

##Make maps to show results
#Polygons with z-scores > 1.96 experience significant positive spatial autocorrelation (SAC)
#Polygons with z-scores < -1.96 experience significant negative spatial autocorrelation (SAC)
#Polygons with z-scores in between -1.96 and 1.96 experience a SAC no different than a polygon exhibiting random SAC
map_LISA_SI <- tm_shape(vriClean) + 
  tm_polygons(col = "Z.Ii_SI", 
              title = "Local Moran's I - Site Index", 
              style = "fixed",
              breaks = c(-Inf, -1.96, 1.96, Inf),
              labels = c("Negative SAC", "Random SAC", "Positive SAC"),
              palette = "RdBu", n = 3, contrast = c(0.19, 0.8),
              midpoint = NA,
              border.alpha = 0.2, colorNA = "black") +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.3))

map_LISA_SI

#----------------------------------------------------------------------------------------------------
###### ELEVATION PPA ######

#if distribution is not random, then the sample isn't a good representation of the study site
#you can critique this in your limitations section
# it probably didn't produce the best interpolated surface as a result, etc. 

##Run PP analysis on point observations
#Choose method and justify (do this for all methods in this study)
#Do NND then KDE to see how it changes over different scales

# tm_shape(vriClean) +
#   tm_polygons(alpha = 0.5) + 
#   tm_shape(elev) + 
#   tm_symbols(size = 0.1, col = "black", alpha = 0.3) + 
#   tm_add_legend(type = "fill", labels = c("grid_code"), col = c(adjustcolor("black", alpha.f = 0.3))) +
#   tm_legend(legend.position = c("left", "bottom")) +
#   tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
#   tm_compass(position=c("left", 0.3))

kma <- elev
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]


#check for and remove duplicated points
#first, finds zero distance among points to see if there are any duplicates
zd <- zerodist(kma)
zd

#if there are duplicates, remove them
kma <- remove.duplicates(kma)

#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(kma)) 

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)


##Average Nearest Neighbour Distance
nearestNeighbour <- nndist(kma.ppp) # for every point we will have a distance from the point to its nearest neighbour

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
n_points = nrow(nearestNeighbour)

nnd = sum(nearestNeighbour$Distance) / n_points

#mean nearest neighbour for random spatial distribution

studyArea <- gArea(vriClean)
pointDensity <- n_points / studyArea

r.nnd = 1 / (2*sqrt(pointDensity))

d.nnd = 1.07453 / sqrt(pointDensity)

c.nnd = 0

R = nnd / r.nnd

SE.NND <- 0.26136 / sqrt(n_points*pointDensity)

z = (nnd-r.nnd) / SE.NND
#IT'S DISPERSED!!!!!!!

#----------------------------------------------------------------------------------------------------
###### CREATE DEM ######

##Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(elev, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(elev)

#IDW Interpolation
#idp = the power value
# you need to change this throughout the code to get the best results
# just play around with it a bit and keep with the same result throughout?
# justify what number you end up choosing?
# idp = 0 means that we don't punish far away values for beng different
# use cross-validation routine to find the best p-value?
# IDW -> the higher the power val, the less significant far away points are in the calculation
# --> Also, it means that the points change very quickly from the origin point. 
P.idw <- gstat::idw(grid_code ~ 1, elev, newdata=grd, idp=4.0)
r       <- raster(P.idw)
r.m     <- mask(r, vriClean)

m3 <-tm_shape(r.m) +
  tm_raster(n=5, palette="Greys", contrast = c(0.2, 1),  
            title="Elevation (m)") +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c("right", "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("right", "top"))

m3 #0

# Leave-one-out validation routine
# this removes one point, and fits an IDW surface to the study area w/out that point
# it does this until every point has been left out once
# this allows us to look at the surface at where that point was and ...
# look at the dif between what we predicted and what the actual was for that left out obs
# by doing this for every point, we can get a sense of the error or influence of leaving one point out
# i.e. the stability of the surface?
# we then have a set of expected and observed points!
IDW.out <- vector(length = length(elev))
for (i in 1:length(elev)) {
  IDW.out[i] <- idw(grid_code ~ 1, elev[-i,], elev[i,], idp=4.0)$var1.pred
}

# Plot the differences
# the red best fit line would be the same aas the black (one-toone) line if the predicted and observed were more similar
# as we make better and better predictions, the red line will slowly tilt to be similar to the black line
# e.g. we predicted 0.036 but we observed 0.030 (see graph)
# as we change the p-value, the predicted and observed should become more similar
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ elev$grid_code, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ elev$grid_code), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
# this spits out a root mean squared error in the units of water variable we are looking at
sqrt( sum((IDW.out - elev$grid_code)^2) / length(elev))


# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(grid_code~1, elev, newdata=grd, idp=4.0)
n   <- length(elev)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(grid_code~1, elev[-i,], newdata=grd, idp=4.0)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Southern California
r_CI <- raster(img.sig, layer="v")
r.m_CI <- mask(r_CI, vriClean)

# Plot the map
# at narrow CIs we are more certain about our estimate
# at wider CIs we are less certain about our estimate
ci_IDW <- tm_shape(r.m_CI) +
  tm_raster(n=3, palette ="Reds", contrast = c(0.2, 1),
            title="95% CI (m)") +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c("right", "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("right", "top"))

tmap_arrange(m3, ci_IDW)

# ##Didn't end up using Kriging since there was an underlying spatial trend in the elevation data: DISPERSION!!
# ##Spatial Interpolation with Kriging
# # Kriging is built on the semivariogram and semivariance!
# # Define the trend model
# f.0 <- as.formula(grid_code ~ 1)
# #Create variogram
# var.smpl <- variogram(f.0, elev, cloud = FALSE) 
# dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
#                           vgm(psill=52000, model="Sph", range=14000, nugget=0))
# 
# # Type ?vgm to see the different model types that you can use in your semivariogram
# # change range, nugget, etc. to get the best fit
# # be able to justify all of this
# plot(var.smpl, dat.fit)
# 
# 
# # Perform the krige interpolation (note the use of the variogram model
# # created in the earlier step)
# dat.krg <- krige(f.0, elev, grd, dat.fit)
# 
# # Convert kriged surface to a raster object for clipping
# r_krig <- raster(dat.krg)
# r.m_krig <- mask(r_krig, vriClean)
# 
# # Plot the map
# elev_DEM <- tm_shape(r.m_krig) +
#   tm_raster(n=5, palette="Greys", contrast = c(0.2, 1),  
#             title="Elevation (m)") +
#   tm_legend(legend.position = c("left", "bottom")) +
#   tm_scale_bar(position=c("right", "bottom"), lwd = 2, text.size=.8) +
#   tm_compass(position=c("right", "top"))
# 
# #This will look at the variance for our kriging surface
# # r_VAR   <- raster(dat.krg, layer="var1.var")
# # r.m_VAR <- mask(r_VAR, vriClean)
# # 
# # tm_shape(r.m_VAR) +
# #   tm_raster(n=3, palette ="Reds", contrast = c(0.2, 1),
# #             title="GVWSA DEM Variance (m)") +
# #   tm_shape(elev) + tm_dots(size=0.2) +
# #   tm_legend(legend.position = c("left", "bottom")) +
# #   tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
# #   tm_compass(position=c("left", 0.3))
# 
# #We can convert our variance into a confidence interval
# r_CI   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
# r.m_CI <- mask(r_CI, vriClean)
# 
# elev_CI <- tm_shape(r.m_CI) + 
#   tm_raster(n=3, palette ="Reds", contrast = c(0.2, 1),
#             title="95% CI (m)") +
#   tm_legend(legend.position = c("left", "bottom")) +
#   tm_scale_bar(position=c("right", "bottom"), lwd = 2, text.size=.8) +
#   tm_compass(position=c("right", "top"))
# 
# tmap_arrange(elev_DEM, elev_CI)

#----------------------------------------------------------------------------------------------------
###### COMBINE VRI DATA WITH ELEVATION ######

##Extract mean value of interpolated surface within each polygon

#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
#Convert your interpolation into a raster and map it:
#r <- raster(r.m_krig)
#sufaceMap <- tm_shape(r.m_krig) + 
#  tm_raster(n=5,palette = "viridis",
#            title="Elev (m)") +
#  tm_shape(elev) + tm_dots(size=0.2)
#sufaceMap

#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average elev for each polygon
vriClean$Elev <- extract(r.m, vriClean, fun = mean)[,1]

#----------------------------------------------------------------------------------------------------
###### LINEAR REGRESSION ######

#Let's say your dataset with both Elev and Height are stored in a dataset called VRI.
#Plot Height and Elev from the VRI dataset you created
plot(vriClean$Site_Index ~ vriClean$Elev)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
VRI.no0 <-  vriClean[which(vriClean$Site_Index > 0), ]
VRI.no0 <-  VRI.no0[which(VRI.no0$Elev > 0), ]

#Now plot the data again
plot(VRI.no0$Site_Index ~ VRI.no0$Elev)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(VRI.no0$Site_Index ~ VRI.no0$Elev)

#Add the regression model to the plot you created
plot(VRI.no0$Site_Index ~ VRI.no0$Elev,
     xlab = "Elevation (m)", ylab = "Site Index (m)", main = "OLS Regression", col='grey30', cex = .8)
abline(lm.model, col = "red", lwd = 2)
abline(h=mean(VRI.no0$Site_Index), col = "blue", lwd = 2)

#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
VRI.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
VRI.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(VRI.no0@data)

#Now, create choropleth map of residuals
map_resid <- tm_shape(VRI.no0) +
  tm_polygons(col = "residuals",
              title = "OLS Residuals (m)",
              style = "fisher",
              palette = "RdBu", n = 5, midpoint = NA,
              contrast = c(0, 0.84), border.alpha = 0.2) +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.3))
map_resid

#----------------------------------------------------------------------------------------------------
###### OLS RESIDUAL SAC ######

### TEST IF RESIDUALS ARE SPATIALLY AUTOCORRELATED OR NOT
##Test assumption of linear regression model and look for independence in the residuals using Global Moran's I
#If there is SAC, then it violates an assumption of linear regression model --> therefore, do GWR

##Defining neighbourhood using using Queen's case
vri.nb2 <- poly2nb(VRI.no0, queen = TRUE)
#Define a network grid
vri.net2 <- nb2lines(vri.nb2, coords=coordinates(VRI.no0))
#Apply a coordinate system
crs(vri.net2) <- crs(VRI.no0)

##Create weights matrix
#If it's a neighbour, it will be 1, if not 0
#For example, if there are two neighbours, the weight of 1 will be distributed ovdfer those two - resulting in a weight of 0.5 for each neighbour
#the default weight style is 'w'?
#Queen's case
vri.lw2 <- nb2listw(vri.nb2, zero.policy = TRUE, style = "W")
print.listw(vri.lw2, zero.policy = TRUE)

##Global Moran's i
#Compute Moran's I Test for site index
#Uses the data from weights matrix, polygon values, and zero.policy settings
#Remember, setting zero.policy = True will bypass the calculation of the polygons with no neighbours
#Setting zero.policy = False will result in a crashed program if there are any polygons that have no neighbours in the data (which there are in this case)
mi_SI <- moran.test(VRI.no0$residuals, vri.lw2, zero.policy = TRUE)
mi_SI

#This function will take ~20 minutes to run
#This calculates the expected range of the spatial autocorrelation given connectivity of the polygons (derived fro the weights matrix)
#moran.range <- function(lw) {
#  wmat <- listw2mat(lw)
#  return(range(eigen((wmat + t(wmat))/2)$values))
#}
#moran.range(vri.lw2)


#This extracts the values needed to calculate the z-score
#Global Moran's I
# FOR SOME REASON THIS DOES NOT WORK --> mi_SI <- mi_SI$estimate[[1]]
mi_SI <- 0.4422919
#Expected Moran's I for a random distribution
#eI_SI = mi_SI$estimate[[2]]
eI_SI <- -0.0001913143
#Variance of values in the dataset
#var_SI <- mi_SI$estimate[[3]]
var_SI <- 0.00007699584
#Calculate Z-score
z_SI <- (mi_SI-eI_SI)/var_SI**0.5

#SIGNIFICANT POSITIVE SAC EXISTS!!!! 95% confident

##Make table to show results
#inlcude other SAC results too?

##Local Moran's i

#this will get a local moran's I value for each polygon
lisa.test_SI <- localmoran(VRI.no0$residuals, vri.lw2)
#extracting the results from the lisa test
VRI.no0$Z.Ii_OLS <- lisa.test_SI[,4]

##Make maps to show results
#Polygons with z-scores > 1.96 experience significant positive spatial autocorrelation (SAC)
#Polygons with z-scores < -1.96 experience significant negative spatial autocorrelation (SAC)
#Polygons with z-scores in between -1.96 and 1.96 experience a SAC no different than a polygon exhibiting random SAC
map_LISA_SI <- tm_shape(VRI.no0) + 
  tm_polygons(col = "Z.Ii_OLS", 
              title = "SAC in OLS Residuals", 
              style = "fixed",
              breaks = c(-Inf, -1.96, 1.96, Inf),
              labels = c("Negative SAC", "Random SAC", "Positive SAC"),
              palette = "RdBu", n = 3, contrast = c(0.19, 0.8),
              midpoint = NA, border.alpha = 0.2, colorNA = "black") +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.3))

map_LISA_SI

#----------------------------------------------------------------------------------------------------
###### GWR ######

#Check assumption of GWR as well

####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
VRI.no0.coords <- sp::coordinates(VRI.no0)
#Observe the result:
head(VRI.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
VRI.no0$X <- VRI.no0.coords[,1]
VRI.no0$Y <- VRI.no0.coords[,2]

###Determine the optimal bandwidth for GWR: this will take a while
#this uses a drop-1 cross validation method
GWRbandwidth <- gwr.sel(VRI.no0$Site_Index ~ VRI.no0$Elev, 
                        data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y),adapt=T)
#--> 0.000628429 (neighbours per unit area)
#--> MOST IMPORTANT PARAMETER IN GWR, IT CONTROLS THE DEGREE OF SMOOTHING IN THE MODEL
#--> FINDS THE OPTIMAL ADAPTIVE NUMBER OF NEIGHBOURS
#----> bandwidth distance changes according to the spatial density of features in the input feature class

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(VRI.no0$Site_Index ~ VRI.no0$Elev, 
                data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model
##Results:
# Call:
#   gwr(formula = VRI.no0$Site_Index ~ VRI.no0$Elev, data = VRI.no0, 
#       coords = cbind(VRI.no0$X, VRI.no0$Y), adapt = GWRbandwidth, 
#       hatmatrix = TRUE, se.fit = TRUE)

# Kernel function: gwr.Gauss 

# Adaptive quantile: 0.0006670305 (about 3 of 5274 data points)

# Summary of GWR coefficient estimates at data points:
#                     Min.     1st Qu.      Median     3rd Qu.        Max.  Global
# X.Intercept. -4.3952e+03  9.0441e+00  2.8891e+01  5.1128e+01  2.6056e+03 26.6797
# VRI.no0.Elev -5.5502e+00 -5.6539e-02 -1.2103e-02  3.0952e-02  6.2844e+00 -0.0066

# Number of data points: 5274 

# Effective number of parameters (residual: 2traceS - traceS'S): 1196.2 
# Effective degrees of freedom (residual: 2traceS - traceS'S): 4077.8 

# Sigma (residual: 2traceS - traceS'S): 4.778823 

# Effective number of parameters (model: traceS): 845.8857 
# Effective degrees of freedom (model: traceS): 4428.114 

# Sigma (model: traceS): 4.585899 --> estimated SD of the residuals. Smaller values is preferable. 
# Sigma (ML): 4.202075 

##AICc --> helpful metric for comparing different regression models
##If a model with AICc value 3 values less than the other model, that means that it is better, compare with OLS!!
#^not an absolute measure of goodness of fit (see R^2 for that)
# AICc (GWR p. 61, eq 2.33; p. 96, eq. 4.21): 32127.68 --> lower values provide better fit to observed data
# AIC (GWR p. 96, eq. 4.22): 30955.33 

# Residual sum of squares: 93125.33  --> the smaller this value, the closer the fir of the GWR model to the observed data
#R2 -> interpreted as the proportion of dependent variable variance accounted for by the regression model
# Quasi-global R2: 0.5759499 


#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
VRI.no0$localr <- results$localR2

#Create choropleth map of r-square values
#Explain why there are negative values
map_r2 <- tm_shape(VRI.no0) +
  tm_polygons(col = "localr",
              title = "GWR r2",
              breaks  = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
              palette = "-viridis", n = 6,
              border.alpha = 0.2) +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.3))
map_r2

#Calculate number of polygons with -R2 values
neg_r2 = length(which(VRI.no0$localr < 0))
tot = length(which(VRI.no0$localr != 0))
neg_r2/tot

vri.r2 <- VRI.no0$localr

#calc num polygons with other r2
num_r2_0.2 = nrow(subset(VRI.no0, localr >= 0 & localr < 0.2))
num_r2_2.4 = nrow(subset(VRI.no0, localr >= .2 & localr < 0.4))
num_r2_2.4/tot
num_r2_4.6 = nrow(subset(VRI.no0, localr >= 0.4 & localr < 0.6))
num_r2_6.8 = nrow(subset(VRI.no0, localr >= 0.6 & localr < 0.8))
num_r2_8.10 = nrow(subset(VRI.no0, localr >= 0.8))

#FIND AREA
VRI.no0$area_sq_m <- area(VRI.no0)
#total area
tot_area = sum(VRI.no0$area_sq_m)
#area with r2 -> 0.2-0.4
area_r2_2.4 = sum(area(subset(VRI.no0, localr >= .2 & localr < 0.4)))
area_r2_2.4/tot_area
#area with negative r2
area_neg_r2 = sum(area(subset(VRI.no0, localr < 0)))
area_neg_r2/tot_area


#Time for more magic. Let's map the coefficients
VRI.no0$coeff <- results$VRI.no0.Elev
#number of stands with negative relationship
num_neg_Coeff = nrow(subset(VRI.no0, coeff < 0))
num_neg_Coeff/tot
#number of stands with positive relationship
num_pos_Coeff = nrow(subset(VRI.no0, coeff > 0))
num_pos_Coeff/tot


#area with positive relationship
area_pos_coeff = sum(area(subset(VRI.no0, coeff > 0)))
area_pos_coeff/tot_area
#area with negative relationship
area_neg_coeff = sum(area(subset(VRI.no0, coeff < 0)))
area_neg_coeff/tot_area

#Create choropleth map of the coefficients
map_coef <- tm_shape(VRI.no0) +
  tm_polygons(col = "coeff",
              title = "GWR Coefficients (m)",
              breaks = c(-4.0, -1.0, 0.0, 1.0, 4.0), midpoint = NA,
              #style = 'fisher', midpoint = NA,
              palette = "BrBG", n = 5,
              border.alpha = 0.2) +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("right", "top"))
map_coef

#----------------------------------------------------------------------------------------------------
###### GWR RESIDUAL SAC ######

##Statistically significant clustering of high and/or low residuals (model under- and overpredictions)
##indicates that the GWR model is misspecified.

#adding predicted and residual column to VRI polygons
VRI.no0$pred <- gwr.model$SDF@data$pred
VRI.no0$gwr.e <- gwr.model$SDF@data$gwr.e

##Defining neighbourhood using using Queen's case
vri.nb2 <- poly2nb(VRI.no0, queen = TRUE)
#Define a network grid
vri.net2 <- nb2lines(vri.nb2, coords=coordinates(VRI.no0))
#Apply a coordinate system
crs(vri.net2) <- crs(VRI.no0)

##Create weights matrix
#If it's a neighbour, it will be 1, if not 0
#For example, if there are two neighbours, the weight of 1 will be distributed ovdfer those two - resulting in a weight of 0.5 for each neighbour
#the default weight style is 'w'?
#Queen's case
vri.lw2 <- nb2listw(vri.nb2, zero.policy = TRUE, style = "W")
print.listw(vri.lw2, zero.policy = TRUE)

##Global Moran's i
#Compute Moran's I Test for site index
#Uses the data from weights matrix, polygon values, and zero.policy settings
#Remember, setting zero.policy = True will bypass the calculation of the polygons with no neighbours
#Setting zero.policy = False will result in a crashed program if there are any polygons that have no neighbours in the data (which there are in this case)
mi_SI <- moran.test(VRI.no0$gwr.e, vri.lw2, zero.policy = TRUE)
mi_SI

#This function will take ~20 minutes to run
#This calculates the expected range of the spatial autocorrelation given connectivity of the polygons (derived fro the weights matrix)
#moran.range <- function(lw) {
#  wmat <- listw2mat(lw)
#  return(range(eigen((wmat + t(wmat))/2)$values))
#}
#moran.range(vri.lw2)

#> moran.range(vri.lw2)
#[1] -1.060660  1.225154

#This extracts the values needed to calculate the z-score
#Global Moran's I
# FOR SOME REASON THIS DOES NOT WORK --> mi_SI <- mi_SI$estimate[[1]]
mi_SI <- -0.02122745
#Expected Moran's I for a random distribution
#eI_SI = mi_SI$estimate[[2]]
eI_SI <- -0.0001913143
#Variance of values in the dataset
#var_SI <- mi_SI$estimate[[3]]
var_SI <- 0.00007698523
#Calculate Z-score
z_SI <- (mi_SI-eI_SI)/var_SI**0.5

#SIGNIFICANT NEGATIVE SAC EXISTS!!!! 95% confident -> Dispersed

##Make table to show results
#inlcude other SAC results too?

##Local Moran's i

#this will get a local moran's I value for each polygon
lisa.test_SI <- localmoran(VRI.no0$gwr.e, vri.lw2)
#extracting the results from the lisa test
VRI.no0$Z.Ii_GWR <- lisa.test_SI[,4]

##Make maps to show results
#Polygons with z-scores > 1.96 experience significant positive spatial autocorrelation (SAC)
#Polygons with z-scores < -1.96 experience significant negative spatial autocorrelation (SAC)
#Polygons with z-scores in between -1.96 and 1.96 experience a SAC no different than a polygon exhibiting random SAC
map_LISA_SI <- tm_shape(VRI.no0) + 
  tm_polygons(col = "Z.Ii_GWR", 
              title = "GWR Residuals", 
              style = "fixed",
              breaks = c(-Inf, -1.96, 1.96, Inf),
              labels = c("Negative SAC", "Random SAC", "Positive SAC"),
              palette = "RdBu", n = 3, contrast = c(0.19, 0.8),
              midpoint = NA,
              border.alpha = 0.2, colorNA = "black") +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_scale_bar(position=c(0.3, "bottom"), lwd = 2, text.size=.8) +
  tm_compass(position=c("left", 0.3))

map_LISA_SI

##plot observed vs predicted to see if GWR was a good prediction
plot(VRI.no0$Site_Index, VRI.no0$pred,
     xlab = "Observed SI (m)", ylab = "Predicted SI (m)", main = "GWR: Observed vs. Predicted SI",
     panel.first=grid())
abline(0, 1, col = "red", lwd = 2)

