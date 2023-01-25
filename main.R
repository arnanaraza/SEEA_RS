### MAIN SCRIPT FOR THE STUDY: 
### Spatial Predictions of Forest Carbon Fluxes for Ecosystem Accounting
### Araza et al. (submitted to journal)

## Preliminaries and global variables
#install.packages('pacman')
pacman::p_load(rgdal,raster,terra,reshape2,ggplot2,Metrics,plyr,dplyr,sf,corrplot,
               rgeos,ggpubr,gridExtra,doParallel,foreach,pdp,ranger,caret,lattice,automap,
               blockCV,vip,hydroGOF,CAST,fastshap,rasterVis,ggnewscale,tidyverse,
               caretEnsemble,viridis,gdalUtilities)

mainDir <- 'C:/SEEA_RS'
outDir <- 'C:/SEEA_RS/results'
dataDir <- 'C:/SEEA_RS/data'

setwd(dataDir)

## Load key functions
source(paste0(mainDir,'/R/Functions.R'))

## Open covariates rasters and reference data table
ref <- read.csv('refLiDAR_deltaAGB_BRA.csv')
covs <- rast('covs_samp.tif')
names(covs) #refer to Table 1 of paper for details

## Mask out non-UNSEEA classes; see Table S3 of paper for reclassification 
## using CCI-Land cover dataset
lc <- covs$lcov10_LC
class1 <- ifel(lc > 49 & lc < 70, 1, NA)#Broadleaved
class2 <- ifel(lc > 69 & lc < 90, 2, NA)#Coniferous
class3 <- ifel(lc > 89 & lc < 110, 3, NA)#Mixed
class4 <- ifel(lc > 159 & lc < 181, 4, NA) #Mangroves
class5 <- ifel(lc > 119 & lc < 131 | 
                 lc > 39 & lc < 50, 5, NA)#Shrubs/Grass/Scubs
lc <- merge(class1,class2,class3,class4,class5)
mgmt <- covs$mgmt_mgmt #Lesiv et al. managed forest data
mgmt <- ifel(mgmt > 30, 6, NA) #Plantations that only overlap classes1-3! 
fmask <- merge(class1,class2)
mgmt <- mask(mgmt,fmask)
lc <- merge(mgmt,lc)
plot(lc)
covs_fmask <- mask(covs,lc)

## Function to graph collinearity among covariates 
Corrplot(covs_fmask) #saves graph to outDir

## Dissimilarity index for HPC runs of AOA
PseudoSamp(ref, covs_fmask) 

## Extract covariates value ----------------------
vt <- ExtractVal(ref, covs_fmask)


#### Model fitting --------------------------
mods <- XVal(vt)


#For validation results per folds
metrics_all <- XvalResults(mods_all)


#### Variogram models -----------------------------------------------
vg <- VG(mods)

