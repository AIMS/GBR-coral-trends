###################################################################
## The following function checks to ensure that all the required ##
## packages are available on the system.                         ##
###################################################################
CoralTrends_checkPackages <- function() {
    require(gdata) # load this first as it masks many
    require(tidyverse)
    require(gtable)
    require(grid)
    require(gridExtra)
    require(xtable)
    library(broom)
    require(rgdal)
    require(rgeos)
    require(sp)
    require(oz)
    require(maps)
    require(mapdata)
    require(ggsn)
    require(scales)
    #require(mapping) # consider replacing this with a self contained function
    require(maptools)
    require(raster)

    require(rstanarm)
    require(coda)
    load(data/primary/whagbr.RData)
    load(data/primary/qld.RData)
}

######################################################################
## The following function is used when issuing a system call from R ##
## (e.g. running xelatex).  It ensures that warnings are captured.  ##
######################################################################
CoralTrends_system <- function (sys.command) {
    ret.val<- system (sys.command, intern=TRUE)

    if (!is.null(attributes(ret.val)$status)) {
        warning (paste('CoralTrends_WARNING', sys.command, " failed | ",paste(ret.val, collapse='\n')))
    }
    ret.val
}


#########################################################################
## The following function appends a log statement into a log file      ##
## parameters:                                                         ##
##     status:    a string indicating either 'FAILURE', 'SUCCESS' or   ##
##                'WARNING'                                            ##
##     logFile:   a character string representation of the log file    ##
##                name (including path relative to current working     ##
##                directory)                                           ##
##     Category:  a character string with a category to appear         ##
##                verbatim in the log                                  ##
##     msg1:      a character string with a message to appear verbatim ##
##                in the log                                           ##
#########################################################################
CoralTrends_log <- function (status, logFile='data/logs/env.log',Category, msg1) {
    options(digits.secs=2)              ## switch to subsecond display
    ## Check if the log file exists, and if it does not, create it
    d=dirname(logFile)
    files <- list.files(d)
    if(!any(grepl(paste0('^',logFile,'$'),files))) system(paste0('touch ',logFile))
    now <- Sys.time()

    msg <- paste0(now, '|',status, ': ', Category, ' ',msg1)
    if( !is.null(msg)){ write(msg,file=paste0(logFile),append=TRUE)}

}

CoralTrends_tryCatch <- function(expr, logFile,Category, expectedClass=NULL, msg=NULL, return=NULL, showWarnings=FALSE) {
    #msg <- paste0(now, '| ', msg)
    max.warnings<-4
    warnings<-0
    W <- NULL
    w.handler <- function(w){ # warning handler
        m<-w$message
        if ((warnings < max.warnings) && (grepl ('CoralTrends_WARNING', m)>0)) {
            CoralTrends_log('WARNING', logFile,Category, paste(warnings, msg, m))
            warnings<<-warnings+1
        }
        invokeRestart("muffleWarning")
    }
    ret <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                    warning = w.handler),warning = W)
    if(!is.atomic(ret$value) && !is.null(ret$value$message)){
        ## An error occurred
        class(ret) <- "try-error"
        CoralTrends_log('WARNING', logFile,Category, paste(msg, ret$value$message))
        if(!is.null(return)) {
            FALSE
        }#else return()
    } else {    #no error check for warning
        CoralTrends_log('INFO', logFile, Category, msg)
        if(!is.null(return)) {
            TRUE
        }
    }
}


CoralTrends_calcPercent = function(x) {
    ifelse(x=='0', 0,
    ifelse(x=='1', 0.05,
    ifelse(x=='1L', 0.025,
    ifelse(x=='1U', 0.075,
    ifelse(x=='2', 0.2,
    ifelse(x=='2L', 0.15,
    ifelse(x=='2U', 0.25,
    ifelse(x=='3', 0.4,
    ifelse(x=='3L', 0.35,
    ifelse(x=='3U', 0.45,
    ifelse(x=='4', 0.625,
    ifelse(x=='4L', 0.5625,
    ifelse(x=='4U', 0.6875,
    ifelse(x=='5', 0.875,
    ifelse(x=='5L',0.8125,0.9375)))))))))))))))
}

COTScategories = function(x) {
    case_when(
        x == 0 ~ 'Zero',
        x > 0 & x < 0.22 ~ 'NO',
        x >= 0.22 & x < 1 ~ 'IO',
        x >= 1 ~ 'AO'
    )
}

CoralTrends_calc3ZoneLocationsOld <- function(x) {
    factor(ifelse(x<= -11.5 & x > -15.4, 'Northern',  #glenn's version is -11.8
           ifelse(x<= -15.4 & x > -20.0, 'Central',
           ifelse(x<= -20.0 & x > -23.92, 'Southern','Outside'))))
}

CoralTrends_calc3ZoneLocations.df <- function(data) {
    data %>% mutate(Zone=factor(ifelse(Latitude<= -11.5 & Latitude > -15.4, 'Northern', #glenn's version is 11.8
                                       ifelse(Latitude<= -15.4 & Latitude > -20.0, 'Central',
                                       ifelse(Latitude<= -20.0 & Latitude > -23.92, 'Southern','Outside')))))
}


CoralTrends_calc3ZoneLocations <- function(manta.sum) {
    load('data/spatial/spatial_3Zone.RData')
    proj4string(spatial_3Zone)=CRS('+proj=longlat +ellps=GRS80 +no_defs')
    ms = manta.sum
    coordinates(ms) = ~Longitude+Latitude
    proj4string(ms)=CRS('+proj=longlat +ellps=GRS80 +no_defs')
    pts=sp:::over(ms, spatial_3Zone)
    manta.sum %>% mutate(zone=pts,Zone=factor(ifelse(zone==1, 'Northern',
                                     ifelse(zone==2, 'Central', 'Southern'))),
                         Location=Zone) %>%
        dplyr::select(-zone)
}

CoralTrends_calc4ZoneLocations <- function(manta.sum.newzones) {
    load('data/spatial/management.RData')
    coordinates(manta.sum.newzones) = ~Longitude+Latitude
    proj4string(manta.sum.newzones)=CRS('+proj=longlat +ellps=GRS80 +no_defs')
    pts=sp:::over(manta.sum.newzones, management)
    pts.wch=if_else(is.na(pts$AREA_DESCR),FALSE,TRUE)
    manta.sum.newzones = cbind(as.data.frame(manta.sum.newzones)[pts.wch,],Zone=pts[pts.wch,'AREA_DESCR'])
    manta.sum.newzones = manta.sum.newzones %>% arrange(desc(Latitude)) %>%
        mutate(Zone=factor(Zone,levels=unique(Zone)))
    save(manta.sum.newzones, file='data/manta.sum.newzones.RData')
    manta.sum.newzones
}


ML_gClip <- function(shp, bb){
    if(class(bb) == "matrix") {
        if (identical(dim(bb), c(2L,2L))) {
            b_poly <- as(raster:::extent(as.vector(t(bb))), "SpatialPolygons")
        } else b_poly = bb
    } else if (class(bb) =='SpatialPolygons') {
        b_poly=bb
    } else b_poly <- as(raster:::extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}


xy2df<-function(xy) {
  data.frame(x=xy$x,y=xy$y)
}
