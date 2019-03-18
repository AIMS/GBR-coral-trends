#########################################################################
## Read in the project specific functions                              ##
## This also includes running a function (CoralTrends_checkPackages()) ##
## that assesses whether all the required packages are present         ##
#########################################################################
source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

#########################################################################
## This codebase includes code for extracting data from AIMS'          ##
## databases.  Currently, this operation can only be performed by      ##
## AIMS staff with adequately credentials. The java applet used to     ##
## perform the data extraction dumps the extraction as a flat csv      ##
## file.  These csv files are available within a publically accessible ##
## data holding hosted by AIMS (link provided *HERE*).                 ##
#########################################################################


#########################################################################
## Read and load the configurations set in parameters/CoralTrends.conf ##
#########################################################################
source('CoralTrends_config.R')

##############################################################
## Generate two alternative zoning configurations           ##
## 1. a three zone configuration based on De'ath et al 2012 ##
##    Zones are stored in:                                  ##
##    - data/spatial/whagbr.RData                           ##
##    - data/spatial/whagbr.n.RData                         ##
##    - data/spatial/whagbr.c.RData                         ##
##    - data/spatial/whagbr.s.RData                         ##
##    - data/spatial/qld.RData                              ##
##############################################################
source('CoralTrends_spatial_3Zone.R')


############################################################################
## Extract and read in the Manta-tow data                                 ##
## and store it under data/primary/manta.csv and data/primary/manta.RData ##
############################################################################
source('CoralTrends_getData_Manta.R')

#################################################
## Process the Manta-tow data                  ##
## - Only include data collected after 1985    ##
## - Convert data into percent cover           ##
## - Summarize the data to the reef/year level ##
## - Assign a zone                             ##
## Data are stored in:                         ##
## - data/processed/manta.sum.RData            ##
## - data/processed/manta.sum.newzones.RData   ##
#################################################
source('CoralTrends_processData_Manta_3Zone.R')

############################
## Generate site maps     ##
## Maps are stored in:    ##
## - figures/MapOfSites.* ##
############################
source('CoralTrends_spatialMap_3Zone.R')

######################################################
## Generate temporal head maps of data availability ##
## Heat maps are stored in:                         ##
## - figures/TemporalHeatMap_3Zones.*               ##
######################################################
source('CoralTrends_temporalHeatMap_3Zones.R')

###############################
## Fit STAN (or INLA) models ##
###############################
source('CoralTrends_stan_manta.R')

##################################
## Disturbance figures and maps ##
##################################
source('CoralTrends_bleaching.R')
source('CoralTrend_cots.R')
source('CoralTrends_cyclones.R')

############################
## Construct trend graphs ##
############################
source('CoralTrend_trend_manta_3Zone.R')

## Disturbance frequency
source('CoralTrends_disturbanceFrequency3Zone.R')

## Map
source('CoralTrends_map_manta.R')

## Disturbance Frequency
source('CoralTrends_disturbanceFrequency_3Zone.R')


## Maps of all the disturbances
source('CoralTrend_spatialFigures.R') 


## Zip up some stuff for Mike
zip('output/NewFigures4Mike.zip', files=c('output/figures/3Zones.pdf',
                                          'output/figures/Disturbances_severe_compilation.pdf',
                                          'output/figures/Fig1mapE_3Zones_new.pdf',
                                          'output/figures/CoralChangeFig.pdf'))



## ----stanAnnualTable, echo=FALSE, results='asis'-------------------------
load(file='output/tables/annualTable.RData')
library(xtable)  
print(xtable(annualTable  %>% as.data.frame,
             caption='Predicted mean ($\\pm$ 95\\% Confidence interval) percent coral cover for selected years (2010--2017) for the GBR, Northern, Central and Southern regions. N indicates the number of observed reefs per location/year combination.'),
      booktabs=TRUE, comment=FALSE,include.rownames=FALSE, caption.placement='top')

## ----stanFullAnnualTable1, echo=FALSE, results='asis'--------------------
load(file='output/tables/fullannualTable.RData')
library(xtable)  
print(xtable(fullannualTable  %>% filter(Location=='Great Barrier Reef') %>% as.data.frame,
             caption='Predicted mean ($\\pm$ 95\\% Confidence interval) percent coral cover for selected years (2010--2017) for the GBR, Northern, Central and Southern regions. N indicates the number of observed reefs per location/year combination.'),
      booktabs=TRUE, comment=FALSE,include.rownames=FALSE, caption.placement='top')

## ----stanFullAnnualTable2, echo=FALSE, results='asis'--------------------
load(file='output/tables/fullannualTable.RData')
library(xtable)  
print(xtable(fullannualTable  %>% filter(Location=='Northern GBR') %>% as.data.frame,
             caption='Predicted mean ($\\pm$ 95\\% Confidence interval) percent coral cover for selected years (2010--2017) for the GBR, Northern, Central and Southern regions. N indicates the number of observed reefs per location/year combination.'),
      booktabs=TRUE, comment=FALSE,include.rownames=FALSE, caption.placement='top')

## ----stanFullAnnualTable3, echo=FALSE, results='asis'--------------------
load(file='output/tables/fullannualTable.RData')
library(xtable)  
print(xtable(fullannualTable  %>% filter(Location=='Central GBR') %>% as.data.frame,
             caption='Predicted mean ($\\pm$ 95\\% Confidence interval) percent coral cover for selected years (2010--2017) for the GBR, Northern, Central and Southern regions. N indicates the number of observed reefs per location/year combination.'),
      booktabs=TRUE, comment=FALSE,include.rownames=FALSE, caption.placement='top')

## ----stanFullAnnualTable4, echo=FALSE, results='asis'--------------------
load(file='output/tables/fullannualTable.RData')
library(xtable)  
print(xtable(fullannualTable  %>% filter(Location=='Southern GBR') %>% as.data.frame,
             caption='Predicted mean ($\\pm$ 95\\% Confidence interval) percent coral cover for selected years (2010--2017) for the GBR, Northern, Central and Southern regions. N indicates the number of observed reefs per location/year combination.'),
      booktabs=TRUE, comment=FALSE,include.rownames=FALSE, caption.placement='top')

## ---- results='markdown'-------------------------------------------------
sessionInfo()

