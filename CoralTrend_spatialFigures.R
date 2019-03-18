source('CoralTrends_functions.R') 
CoralTrends_checkPackages()



library(rgeos)
library(sp)
library(tidyverse)
library(mgcv)
library(INLA)
library(oz)
library(maps)
library(mapdata)
library(mapping)
library(ggsn)
library(scales)
library(grid)
library(broom)
require(maps)
require(maptools)
library(gtable)
source('CoralTrends_functions.R')

ZONES=3

#########################################################################
## This script is used to generate graphical overviews of the spatial  ##
## sampling design, COTS, Cyclones, Bleaching data for the period 2013 ##
## -- 2017.                                                            ##
## There are two possible approaches:                                  ##
## 1. simply plot points reflecting categorical observed data          ##
## 2. underlay the points with a spatial surface.                      ##
#########################################################################


gClip <- function(shp, bb, type='Intersection'){
    require(raster)
    require(rgeos)
    if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
    else b_poly <- as(extent(bb), "SpatialPolygons")
    if(type=='Intersection') return(gIntersection(shp, b_poly, byid = T))
    if(type=='Difference') return(gDifference(b_poly,shp, byid = T))
}

display_Long = function(data, Long=0) {
    data %>% mutate(long=long+Long)
    }



## Approach 1.
load('data/spatial/qld.RData')
load(file='data/spatial/management.RData')
load(file='data/spatial/spatial_3Zone.RData')
if (ZONES==3) management=spatial_3Zone
bb = bbox(management)
management.df = fortify(management)
tops=management.df %>% filter(lat > max(lat)-0.1) %>% summarize(min=min(long),max=max(long), lat=max(lat))

if (ZONES==3) {
    Zones=c('Northern','Central','Southern')
} else {
    Zones=gsub(' Management Area','',rev(management@data$AREA_DESC))
    Zones=gsub('/','\n',Zones)
    Zones=gsub(' ','\n',Zones)
}

## Sampling design
if (ZONES==4) {
    load('data/processed/manta.sum.newzones.RData')
    manta.sum = manta.sum.newzones %>%
                                        #filter(REPORT_YEAR > 2012) %>%
        group_by(P_CODE.mod,Zone,REEF_NAME,REEF_ID) %>%
        summarize(Latitude=mean(Latitude,na.rm=TRUE), Longitude=mean(Longitude,na.rm=TRUE),Tows=mean(Tows,na.rm=TRUE))
} else {
    load('data/processed/manta.sum.RData')
    manta.sum = manta.sum %>%
                                        #filter(REPORT_YEAR > 2012) %>%
        group_by(P_CODE.mod,Zone,REEF_NAME,REEF_ID) %>%
        summarize(Latitude=mean(Latitude,na.rm=TRUE), Longitude=mean(Longitude,na.rm=TRUE),Tows=mean(Tows,na.rm=TRUE))
}


gbr.sp<- readShapeSpatial("~/Work/Resources/GIS/Features/Great_Barrier_Reef_Features.shp",
                          proj4string = CRS('+proj=longlat +ellps=WGS84'),repair=TRUE,force_ring=T,verbose=TRUE)
## Clipped coast
gbr.sp = raster:::crop(gbr.sp, raster::extent(bb))
gbr.fort <- tidy(gbr.sp, region='FEAT_NAME')


grd.y = data.frame(x=c(140,140,140,140,140,140,140),
                 xend=c(146,146,148,150,152,154,154),
                 y=c(-11,-13,-15,-17,-19,-21,-23),
                 yend=c(-11,-13,-15,-17,-19,-21,-23))
grd.x = data.frame(x=c(142,144,146,148,150,152,154),
                 xend=c(142,144,146,148,150,152,154),
                 y=c(-9,-9,-11,-15,-17,-19,-21),
                 yend=c(-30,-30,-30,-30,-30,-30,-30))
g.base = ggplot(manta.sum, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long, group=group),fill='grey', color='grey70') +
    geom_polygon(data=fortify(qld), aes(y=lat, x=long, group=group),fill='grey', color='grey40') +
    geom_segment(data=grd.y, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80')+
    geom_segment(data=grd.x, aes(x=x,y=y,xend=xend,yend=yend), size=0.1, color='grey80')+
    annotate(geom='text', x=-Inf,y=seq(-23,-11,by=2), label=degreeLabelsNS(seq(-23,-11,by=2)), hjust=-0.1,parse=TRUE, size=2) +
    annotate(geom='text', y=-Inf,x=seq(142,154,by=2), label=degreeLabelsEW(seq(142,154,by=2)),vjust=-0.2,parse=TRUE,size=2) +
    geom_point(shape=21,fill='yellow',color='black') +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.5) +
    geom_point(data=layers:::towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long)) +
    geom_text(data=layers:::towns12 %>% filter(town %in% c('Cooktown','Cairns','Townsville','Mackay','Rockhampton')), aes(y=lat, x=long, label=town), hjust=1.1, size=3)
if (ZONES==4) {
    g.base = g.base + annotate(geom='text',y=-11.8,x=144.7,label=Zones[1], angle=-90, hjust=0.5,vjust=-0.5) +
                                        #annotate(geom='segment',y=-12,x=145.5,xend=145,yend=-12) +
        annotate(geom='text',y=-16,x=146,label=Zones[2], angle=-68.5, hjust=0.5,vjust=-0.5) +
                                        #annotate(geom='segment',y=-15.6,x=147,xend=146.5,yend=-15.6) +
        annotate(geom='text',y=-19.1,x=149,label=Zones[3], angle=-31, hjust=0.5,vjust=-0.5) +
        annotate(geom='text',y=-20.4,x=152.0,label='Mackay', angle=-30.8, hjust=0.5,vjust=-0.5) +
        annotate(geom='text',y=-22,x=153.2,label='Capricorn', angle=-73, hjust=0.5,vjust=-0.5)
} else {
    g.base = g.base + annotate(geom='text',y=-13,x=145.2,label=Zones[1], angle=-76, hjust=0.5,vjust=-0.5) +
        annotate(geom='text',y=-18,x=147.8,label=Zones[2], angle=-31, hjust=0.5,vjust=-0.5) +
        annotate(geom='text',y=-22,x=153.2,label=Zones[3], angle=-73, hjust=0.5,vjust=-0.5)
}
g.base=g.base + scale_x_continuous(breaks=seq(142,154,by=2), position = 'bottom') +
    scale_y_continuous(breaks=seq(-25,-10,by=2)) +
    coord_equal(ylim=bbox(management)[2,],xlim=c(bbox(management)[1,],bbox(management)[1,]+20)) +
    theme_minimal() +
    theme(legend.position=c(0.05,0.10), legend.justification=c(0,0),
          legend.background = element_rect(color='black', fill='white'),
          axis.title=element_blank(),
          panel.border = element_rect(color='black',fill=NA),
          axis.text=element_blank(),
                                        #axis.ticks.length = unit(c(-0.5,-0.9),'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g.base


# calculate slope (for the purpose of having labels with text aligned along a slope)
if (ZONES==4) dat=management.df%>% filter(lat< -15 & lat> -16, long>146) #%>% filter(lat==first(lat) | lat==last(lat))
if (ZONES==3) dat=management.df%>% filter(lat< -15 & lat> -16, long>145.5) #%>% filter(lat==first(lat) | lat==last(lat))
(coefs=coef(lm(long~lat,dat)))
strt = c(145.9,-14.5)
ff=function(strt,coefs,length) {
    x=strt[1]
    y=strt[2]
    dat=data.frame(x=x,y=y)
    for (i in 2:length) {
        y=y-0.7
        x=coefs[1] + y*coefs[2]
        dat=rbind(dat, data.frame(x=x[[1]], y=y[[1]]))
    }
    dat
}

## Bleaching
#load(file='data/processed/bleaching.2016.RData')
load(file='data/processed/bleaching.2016_3Zone.RData')
bleaching.2016 = bleaching.2016 %>%
    group_by(REEF_NAME) %>%
    summarize(Latitude=mean(Latitude, na.rm=TRUE), Longitude=mean(Longitude, na.rm=TRUE), BLEACHINGcat=max(as.numeric(as.character(BLEACHINGcat)),na.rm=TRUE)) %>%
    mutate(y=as.numeric(BLEACHINGcat)-1, Cat=BLEACHINGcat, Cat=factor(Cat, levels=0:4, labels=paste('BL',0:4)))
nCat.bleaching = levels(bleaching.2016$Cat)

clrs.bleaching = c('#FFFFFF',brewer_pal(palette='Reds')(5)[-1])
shapes.bleaching = c(21:25)
legd.bleaching <- legendGrob(c("0","1","2","3","4"), pch=21,gp=gpar(col = 'black', fill = clrs.bleaching),vgap=unit(0.15,'lines'))
colLim=levels(bleaching.2016$Cat)
lgnd.dat = ff(strt,coefs,length=5)
lgnd.dat$lab = c('0','1','2','3','4')
lgnd.dat$Cat=factor(paste('BL',lgnd.dat$lab))
g.bleaching = ggplot(manta.sum, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=bleaching.2016 %>% mutate(Longitude=Longitude),aes(fill=Cat,shape=Cat)) +
    #geom_point(data=bleaching.2016 %>% mutate(Longitude=Longitude),aes(size=Cat),fill=clrs.bleaching[5],shape=21,alpha=0.4) +
    scale_fill_manual(limit=colLim,values=clrs.bleaching)+
    scale_shape_manual(limit=colLim,values=shapes.bleaching)+
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Bleaching'), vjust=0) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(data=lgnd.dat, aes(y=y,x=x+0.4, fill=Cat,shape=Cat), size=3) +
    #geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, size=Cat), shape=21, fill=clrs.bleaching[5]) +
    geom_text(data=lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

g.bleaching=ggplotGrob(g.bleaching)
g.bleaching.grob.panel=g.bleaching[[1]][[6]]
#legd.bleaching <- legendGrob(c("0","1","2","3","4"), pch=21,
#                             gp=gpar(col = 'black', fill = clrs[1:5]),vgap=unit(0.15,'lines'))
                                        #annotation_custom(grob=legd.bleaching,xmin=150,xmax=151,ymax=-11,ymin=-13) +
shiftacross = 7#6
g = g.base + annotation_custom(grob=g.bleaching.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf) 


## COTS
cots <- read.csv('data/primary/cots.csv',strip.white=TRUE)
#load('data/manta.sum.newzones.RData')
#cots.2016 = cots %>% filter(REEF_NAME %in% as.character(unique(manta.sum.newzones$REEF_NAME))) %>%
cots.2016 = cots %>% filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) %>%
    mutate(COTScat = COTScategories(AVGOFMEAN_COTS), COTScat = factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    #left_join(manta.sum.newzones %>% dplyr:::select(REEF_NAME,Zone) %>% distinct) %>%
    left_join(manta.sum %>% dplyr:::select(REEF_NAME,Zone) %>% distinct) %>%
    dplyr:::select(REEF_NAME,COTScat, Latitude=REEF_LAT, Longitude=REEF_LONG, REPORT_YEAR, Zone) %>%
    filter(REPORT_YEAR>2012) %>%
    mutate(y=as.numeric(COTScat)-1, COTScat=factor(y, levels=0:3,labels=c('Zero','NO','IO','AO')),
           Cat=as.numeric(as.character(factor(COTScat, labels=0:3)))) %>%
    group_by(REEF_NAME) %>%
    summarize(Latitude=mean(Latitude, na.rm=TRUE), Longitude=mean(Longitude, na.rm=TRUE), Cat=max(Cat,na.rm=TRUE))%>%
    mutate(Cat=factor(Cat, levels=0:3, labels=paste('CO',c('Zero','NO','IO','AO')))) %>%
    arrange(Cat)

clrs.cots = c('#FFFFFF',brewer_pal(palette='Greens')(4)[-1])
shapes.cots = c(21:24)
legd.cots <- legendGrob(c("0","1","2","3"), pch=21,gp=gpar(col = 'black', fill = clrs.cots),vgap=unit(0.15,'lines'))
colLim=levels(cots.2016$Cat)
lgnd.dat = ff(strt,coefs,length=4)
lgnd.dat$lab = c('Zero','NO','IO','AO')
lgnd.dat$Cat=factor(paste('CO',lgnd.dat$lab))
g.cots = ggplot(manta.sum, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=cots.2016 %>% mutate(Longitude=Longitude),aes(fill=Cat,shape=Cat)) +
    scale_fill_manual(limit=colLim,values=clrs.cots)+
    scale_shape_manual(limit=colLim,values=shapes.cots)+
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='COTS'), vjust=0) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(data=lgnd.dat, aes(y=y,x=x+0.4, fill=Cat,shape=Cat),size=3) +
    geom_text(data=lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  

g.cots=ggplotGrob(g.cots)
g.cots.grob.panel=g.cots[[1]][[6]]
shiftacross = shiftacross+7 #11
g = g + annotation_custom(grob=g.cots.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf)

## Cyclones
cyclones = read.csv('data/primary/170909 Cyclone wave data from Marji.csv', strip.white=TRUE)
cyclones.2016 = cyclones %>%
    dplyr:::rename(REEF_NAME=Reef,A_SECTOR=TMP_sector,SHELF=Shelf,REEF_LAT=lat, REEF_LONG=long_) %>%
    filter(!REEF_NAME=='') %>%                                                   #remove cases with missing reef name
    ##dplyr:::select(-FID_,-max) %>%                                               #remove the FID_ field which has no data and max (which I am not sure what to do with)
    dplyr:::select(-Project,-gbrmpa_sector,-full_reef_id,-gazetted_name,-GBRMPA_ID) %>% #remote new and unnecessary fields
    ##gather(key=REPORT_YEAR, value=S, -A_SECTOR:-SAMPLE_CLA) %>%                  #melt the data by year columns
    gather(key=REPORT_YEAR, value=S, -A_SECTOR:-REEF_LONG) %>%                    #melt the data by year columns
    mutate(REPORT_YEAR=as.numeric(as.character(gsub('X','',REPORT_YEAR)))) %>%   #generate a REPORT_YEAR field that is a numeric year
    filter(REPORT_YEAR>2012) %>%
    mutate(CYCLONEcat = S, Cat=CYCLONEcat) %>%                                   #generate a CYCLONEcat field that has the cyclone category
    group_by(REEF_NAME) %>%
    summarize(Latitude=mean(REEF_LAT, na.rm=TRUE), Longitude=mean(REEF_LONG, na.rm=TRUE),Cat=max(Cat, na.rm=TRUE)) %>%
    mutate(Cat=factor(Cat, levels=0:3, labels=paste('CY',c("0","1","2","3"))))

clrs.cyclones = c('#FFFFFF',brewer_pal(palette='Blues')(4)[-1])
shapes.cyclones = c(21:24)
legd.cyclones <- legendGrob(c("0","1","2","3"), pch=21,gp=gpar(col = 'black', fill = clrs.cyclones),vgap=unit(0.15,'lines'))
colLim=levels(cyclones.2016$Cat)
lgnd.dat = ff(strt,coefs,length=4)
lgnd.dat$lab = c("Hs<4*m","Hs>=4*m~'&'<6*m","Hs>=6*m~'&'<8*m", "Hs>=8*m")
lgnd.dat$Cat=factor(paste('CY',0:3))
g.cyclones = ggplot(manta.sum, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=cyclones.2016 %>% mutate(Longitude=Longitude),aes(fill=Cat,shape=Cat)) +
    scale_fill_manual(limit=colLim,values=clrs.cyclones)+
    scale_shape_manual(limit=colLim,values=shapes.cyclones)+
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Cyclones'), vjust=0) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(data=lgnd.dat, aes(y=y,x=x+0.4, fill=Cat, shape=Cat),size=3) +
    geom_text(data=lgnd.dat, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0, parse=TRUE) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())   
g.cyclones=ggplotGrob(g.cyclones)
g.cyclones.grob.panel=g.cyclones[[1]][[6]]
shiftacross = shiftacross + 7 #16
g = g + annotation_custom(grob=g.cyclones.grob.panel, xmin=bb[1,1]+shiftacross, xmax=bb[1,2]+shiftacross, ymin=-Inf,ymax=Inf)


## Coral change - this used to be on the right hand side of the top panel as an average annual change since 2012.
## Mike would now like a separate panel that illustrates the annual change in coral cover for each of the years 2012 - 2018
## so as to illustrate annual change each year

## ## Mike supplied the following - Need to discuss this with him
## ann.change<-read.csv(file='parameters/post2012.csv',strip.white=T)
## ann.change.ave<-ann.change %>% group_by(GBRMPA,REGION,A_SECTOR,SHELF,REEF_ID,REEF_NAME) %>% dplyr:::summarize(mean=mean(annual.change))
## ann.change.ave$REGION<-factor(ann.change.ave$REGION,levels=c("Northern","Central","Southern"))
## ann.change.ave$GBRMPA<-factor(ann.change.ave$GBRMPA,levels=c("Far Northern","Cairns/Cooktown","Townsville/Whitsunday","Mackay/Capricorn"))
## ann.change.se<-ann.change %>% group_by(GBRMPA,REGION,A_SECTOR,SHELF,REEF_ID,REEF_NAME) %>% dplyr:::summarize(SE=sd(annual.change)/sqrt(length(annual.change)))
## ann.change.se$REGION<-factor(ann.change.ave$REGION,levels=c("Northern","Central","Southern"))
## ann.change.se$GBRMPA<-factor(ann.change.ave$GBRMPA,levels=c("Far Northern","Cairns/Cooktown","Townsville/Whitsunday","Mackay/Capricorn"))
## ann.change.ave$SE<-ann.change.se$SE
## ann.change.ave$lower<-ann.change.ave$mean-ann.change.ave$SE
## ann.change.ave$upper<-ann.change.ave$mean+ann.change.ave$SE
## coral.2016 = ann.change.ave %>% left_join(manta.sum %>% dplyr:::select(REEF_NAME,Latitude,Longitude) %>% distinct)

## OR we could recreate using methodology that is more internally consistent.
load('data/processed/manta.sum.RData')
manta.change=manta.sum %>% filter(REPORT_YEAR>=2012) %>%
    group_by(Location,Zone,REEF_NAME,Longitude,Latitude) %>%
    mutate(Cover.lag = lag(Cover),
           Change = 100*(Cover-Cover.lag),#/Cover.lag,
           Time.lag=lag(REPORT_YEAR),
           Time.change=REPORT_YEAR-Time.lag,
           Change.annual=Change/Time.change) %>%
    filter(REPORT_YEAR>2012) %>%
    group_by(Location,Zone,REEF_NAME,Longitude,Latitude) %>%
    summarize(mean=mean(Change.annual, na.rm=TRUE)) %>%
    ungroup %>%
    filter(!is.na(mean))
coral.2016 = manta.change

ann.change = manta.sum %>% filter(REPORT_YEAR==2012) %>% droplevels %>%
    full_join(manta.sum %>% group_by(REEF_NAME) %>% filter(REPORT_YEAR==max(REPORT_YEAR)) %>% dplyr::select(REEF_NAME, Cover2=Cover, REPORT_YEAR2=REPORT_YEAR)) %>%
    ungroup %>% droplevels %>%
    mutate(change=100*(Cover2-Cover),
           number.of.years=REPORT_YEAR2-REPORT_YEAR,
           annual.change=change/number.of.years)

clrs.coral = c('#FF0000','#FF00FF')
legd.coral <- legendGrob(c("Decrease","Increase"), pch=21,gp=gpar(col = 'black', fill = clrs.coral),vgap=unit(0.15,'lines'))
#colLim=levels(coral.2016$Cat)
lgnd.dat = data.frame(x=c(149,149),y=c(-16,-16.7)) #ff(strt,coefs,length=2)
lgnd.dat$lab = c("Decrease","Increase")
lgnd.dat$Cat=factor(c(FALSE,TRUE))
lgnd.dat.size = ff(lgnd.dat[2,1:2],coefs,length=5)
lgnd.dat.size = ff(strt,coefs,length=5)
lgnd.dat.size$lab = format(seq(2.5,12.5,by=2.5),digits=2)
lgnd.dat.size$size = seq(2.5,12.5,by=2.5)

g.coral = ggplot(manta.sum, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=coral.2016 %>% mutate(Longitude=Longitude),aes(fill=mean>0, size=abs(mean)),shape=21,alpha=0.4) +
    scale_fill_manual('',breaks=c(FALSE,TRUE), values=c('red','green'), labels=c('<0','>0')) +
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.1, x=long,label='Coral change'), vjust=0) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
    geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=TRUE) +
    geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
    geom_point(data=lgnd.dat.size, aes(y=y,x=x+0.4, size=size), shape=21) +
    geom_text(data=lgnd.dat.size, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0, parse=TRUE) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())    +
    scale_radius()
g.coral

g.coral=ggplotGrob(g.coral)
g.coral.grob.panel=g.coral[[1]][[6]]
#g = g + annotation_custom(grob=g.coral.grob.panel, xmin=bb[1,1]+21, xmax=bb[1,2]+21, ymin=-Inf,ymax=Inf)

## Add the Australia inset
coord.df <- map_data("world2", "Australia", exact=FALSE,
                  xlim=c(110,155),ylim=c(-45,-5),
                  boundary=FALSE,
                  interior=TRUE, as.polygon=TRUE)
coord <- maps:::map.poly("worldHires", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=TRUE, as.polygon=TRUE)
IDs <- sapply(strsplit(coord$names, ":"), function(x) x[1])
coord.sp <- map2SpatialPolygons(coord,IDs=IDs)
bb.4=bb
bb.4[1,1] = bb[1,1]-(bb[1,2] - bb[1,1])*0.1
bb.4[2,1] = bb[2,1]-(bb[2,2] - bb[2,1])*0.1
coord.sp = gClip(coord.sp, bb.4)

aus=ggplot(coord.df, aes(y=lat,x=long,group=group)) + geom_polygon(fill='white',color='black') + coord_map() +
    geom_polygon(data=fortify(coord.sp), aes(y=lat, x=long, group=group),fill='grey',color='black', size=0.2) +
    geom_polygon(data=fortify(management), aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    #annotate(geom='rect', xmin=bb[1,1],xmax=bb[1,2], ymin=bb[2,1],ymax=bb[2,2], fill=NA, color='black') +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_rect(color='black'))
aus.grob=ggplotGrob(aus)
g = g + annotation_custom(grob=aus.grob,xmax=Inf,ymax=Inf,xmin=170,ymin=-15) 


## Add the scalebar and north arrow
gc.1=c(142.5,-24)
gc.2=gcDestination(lon=gc.1[1], lat=gc.1[2],bearing=90, dist=500, model='WGS84', Vincenty=FALSE)
g = g + ggsn:::scalebar(x.min=gc.1[1],y.min=gc.1[2],x.max=gc.2[1],y.max=gc.1[2]+0.5,dist=250, dd2km=TRUE, model='WGS84',st.size=2, height=0.5,st.dist=0.5) +
    ggsn:::north(x.min=mean(c(gc.1[1],gc.2[1]))-0.5, x.max=mean(c(gc.1[1],gc.2[1]))+0.5, y.min=gc.1[2]+0.5,y.max=gc.1[2]+2, scale=1) 
g1=g




## Cyclone paths
##Cyclone paths
dylan.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/dylan_active_line.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
dylan.sp <- spTransform(dylan.sp, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
dylan.fort <- tidy(dylan.sp)
dylan.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/poly_4mw_dylan.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
dylan.4m <- spTransform(dylan.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
dylan.4m.fort <- tidy(dylan.4m)
dylan.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w6m_dylan.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
dylan.6m <- spTransform(dylan.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
dylan.6m.fort <- tidy(dylan.6m)
#------------------------


ita.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/ita_active_line.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
ita.sp <- spTransform(ita.sp, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ita.fort <- tidy(ita.sp)
ita.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/poly_4mw_ita14_use.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
ita.4m <- spTransform(ita.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ita.4m.fort <- tidy(ita.4m)
ita.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w6m_ita.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
ita.6m <- spTransform(ita.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ita.6m.fort <- tidy(ita.6m)
ita.8m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w8m_ita.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
ita.8m <- spTransform(ita.8m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ita.8m.fort <- tidy(ita.8m)

#--------
oswald.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/oswald_active_line.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
oswald.fort <- tidy(oswald.sp)
oswald.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w4m_oswald.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
oswald.4m <- spTransform(oswald.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
oswald.4m.fort <- tidy(oswald.4m)
oswald.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w6m_oswald.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
oswald.6m <- spTransform(oswald.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
oswald.6m.fort <- tidy(oswald.6m)
                                        #----------------------

marcia.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/marcia_line.shp",
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
marcia.fort <- tidy(marcia.sp)
marcia.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/poly_4mw_marcia.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
marcia.4m <- spTransform(marcia.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
marcia.4m.fort <- tidy(marcia.4m)
marcia.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w6m_marcia.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
marcia.6m <- spTransform(marcia.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
marcia.6m.fort <- tidy(marcia.6m)
marcia.8m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w8m_marcia.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
marcia.8m <- spTransform(marcia.8m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
marcia.8m.fort <- tidy(marcia.8m)

#--------------------
nathan.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/nathan_line.shp",
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
nathan.fort <- tidy(nathan.sp)
nathan.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/poly_4mw_nathan.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
nathan.4m <- spTransform(nathan.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
nathan.4m.fort <- tidy(nathan.4m)
nathan.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w6m_n15.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
nathan.6m <- spTransform(nathan.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
nathan.6m.fort <- tidy(nathan.6m)
nathan.8m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/gbr_w8m_n15.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
nathan.8m <- spTransform(nathan.8m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
nathan.8m.fort <- tidy(nathan.8m)

#--------------------
debbie.sp<- readShapeSpatial("parameters/GIS/cyclone_tracks_for_figure/tc_debbie_2017_track_line.shp",
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
debbie.fort <- tidy(debbie.sp)
debbie.4m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/w4m_envelope_debbie17.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
debbie.4m <- spTransform(debbie.4m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
debbie.4m.fort <- tidy(debbie.4m)
debbie.6m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/w6m_envelope_debbie17.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
debbie.6m <- spTransform(debbie.6m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
debbie.6m.fort <- tidy(debbie.6m)
debbie.8m<- readShapeSpatial("parameters/GIS/cyclone_wave_envelopes/w8m_envelope_debbie17.shp",
                            proj4string = CRS('+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
debbie.8m <- spTransform(debbie.8m, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
debbie.8m.fort <- tidy(debbie.8m)



## g2=ggplot(data=NULL) +
##     geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long, group=group),fill='grey', color='grey70') +
##     #geom_polygon(data=fortify(qld), aes(y=lat, x=long, group=group),fill='grey', color='grey40') +
##     geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
##     geom_polygon(data=oswald.4m.fort, aes(y=lat, x=long+0, group=group, alpha='Hs 4'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=oswald.6m.fort, aes(y=lat, x=long+0, group=group, alpha='Hs 6'), fill='blue',show.legend = FALSE) +
##     geom_path(data=oswald.fort, aes(y=lat, x=long), color='blue') +
##     geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Oswald (Jan 2013)'), vjust=0) +

##     geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+10, group=group),fill='grey', color='grey70') +
##     geom_polygon(data=management.df, aes(y=lat, x=long+10, group=group),fill=NA,color='black', size=0.2) +  
##     geom_polygon(data=dylan.4m.fort, aes(y=lat, x=long+10, group=group, alpha='Hs 4'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=dylan.6m.fort, aes(y=lat, x=long+10, group=group, alpha='Hs 6'), fill='blue',show.legend = FALSE) +
##     geom_path(data=dylan.fort, aes(y=lat, x=long+10), color='blue') +
##     geom_text(data=tops %>% mutate(min=min+10, max=max+10, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Dylan (Jan 2014)'), vjust=0) +

##     geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+20, group=group),fill='grey', color='grey70') +
##     geom_polygon(data=management.df, aes(y=lat, x=long+20, group=group),fill=NA,color='black', size=0.2) +      
##     geom_polygon(data=ita.4m.fort, aes(y=lat, x=long+20, group=group, alpha='Hs 4'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=ita.6m.fort, aes(y=lat, x=long+20, group=group, alpha='Hs 6'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=ita.8m.fort, aes(y=lat, x=long+20, group=group, alpha='Hs 8'), fill='blue',show.legend = FALSE) +
##     geom_path(data=ita.fort %>% filter(long<150 | lat < -17), aes(y=lat, x=long+20), color='blue') +
##     geom_text(data=tops %>% mutate(min=min+20, max=max+20, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Ita (Apr 2014)'), vjust=0) +

##     geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+30, group=group),fill='grey', color='grey70') +    
##     geom_polygon(data=management.df, aes(y=lat, x=long+30, group=group),fill=NA,color='black', size=0.2) +  
##     geom_path(data=marcia.fort %>% filter(long<155), aes(y=lat, x=long+30), color='blue') +
##     geom_polygon(data=marcia.4m.fort, aes(y=lat, x=long+30, group=group, alpha='Hs 4'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=marcia.6m.fort, aes(y=lat, x=long+30, group=group, alpha='Hs 6'), fill='blue',show.legend = FALSE) +
##     geom_polygon(data=marcia.8m.fort, aes(y=lat, x=long+30, group=group, alpha='Hs 8'), fill='blue',show.legend = FALSE) +
##     geom_text(data=tops %>% mutate(min=min+30, max=max+30, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Marcia (Feb 2015)'), vjust=0) +
    
##     geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+40, group=group),fill='grey', color='grey70') +
##     geom_polygon(data=management.df, aes(y=lat, x=long+40, group=group),fill=NA,color='black', size=0.2) +  
##     geom_path(data=nathan.fort, aes(y=lat, x=long+40), color='blue') +
##     geom_polygon(data=nathan.4m.fort, aes(y=lat, x=long+40, group=group, alpha='Hs 4'), fill='blue') +
##     geom_polygon(data=nathan.6m.fort, aes(y=lat, x=long+40, group=group, alpha='Hs 6'), fill='blue') +
##     geom_polygon(data=nathan.8m.fort, aes(y=lat, x=long+40, group=group, alpha='Hs 8'), fill='blue') +
##     geom_text(data=tops %>% mutate(min=min+40, max=max+40, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Nathan (Mar 2015)'), vjust=0) +

##     scale_alpha_manual('Wave exposure', breaks = c('Hs 0','Hs 4', 'Hs 6', 'Hs 8'), labels=c(expression(Hs<4*m),expression(Hs>=4*m~and<6*m), expression(Hs>=6*m~and<8*m), expression(Hs>=8*m)), values=c(0,0.1,0.2,0.3), limits=c('Hs 0','Hs 4', 'Hs 6', 'Hs 8'))+
##     scale_x_continuous(breaks=seq(142,154,by=2), position = 'top') +
##     scale_y_continuous(breaks=seq(-25,-10,by=2)) +
##     coord_equal(ylim=bbox(management)[2,]+c(0,1),xlim=c(bbox(management)[1,],bbox(management)[1,]+40)) +
##     theme_minimal() +
##     theme(legend.position=c(0.01,0.01), legend.justification=c(0,0),
##           legend.background = element_blank(),
##           legend.text.align = 0,legend.key=element_rect(color='black',fill=NA),
##           axis.title.y=element_blank(), axis.text.y=element_blank(),
##           axis.title.x=element_blank(), axis.text.x=element_blank(),
##           panel.border = element_rect(color='black',fill=NA),
##           panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#legd.cyclones <- legendGrob(c(expression(Hs<4*m),expression(Hs>=4*m~'&'<6*m), expression(Hs>=6*m~'&'<8*m), expression(Hs>=8*m)), pch=21,
#                             gp=gpar(col = 'black', fill = clrs.cyclones),vgap=unit(0.15,'lines'))

reef_offset = 8
g2=ggplot(data=NULL) +
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long, group=group),fill='grey', color='grey70') +
    #geom_polygon(data=fortify(qld), aes(y=lat, x=long, group=group),fill='grey', color='grey40') +
    geom_polygon(data=oswald.4m.fort, aes(y=lat, x=long+0, group=group), fill=clrs.cyclones[2],show.legend = FALSE) +
    geom_polygon(data=oswald.6m.fort, aes(y=lat, x=long+0, group=group), fill=clrs.cyclones[3],show.legend = FALSE) +
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Oswald (Jan 2013)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_path(data=oswald.fort, aes(y=lat, x=long), color='blue') +
    
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+reef_offset, group=group),fill='grey', color='grey70') +
    geom_polygon(data=dylan.4m.fort, aes(y=lat, x=long+reef_offset, group=group), fill=clrs.cyclones[2],show.legend = FALSE) +
    geom_polygon(data=dylan.6m.fort, aes(y=lat, x=long+reef_offset, group=group), fill=clrs.cyclones[3],show.legend = FALSE) +
    geom_text(data=tops %>% mutate(min=min+reef_offset, max=max+reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Dylan (Jan 2014)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+reef_offset, group=group),fill=NA,color='black', size=0.2) +  
    geom_path(data=dylan.fort %>% filter(long<151.5), aes(y=lat, x=long+reef_offset), color='blue') +
    
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+2*reef_offset, group=group),fill='grey', color='grey70') +
    geom_polygon(data=ita.4m.fort, aes(y=lat, x=long+2*reef_offset, group=group), fill=clrs.cyclones[2],show.legend = FALSE) +
    geom_polygon(data=ita.6m.fort, aes(y=lat, x=long+2*reef_offset, group=group), fill=clrs.cyclones[3],show.legend = FALSE) +
    geom_polygon(data=ita.8m.fort, aes(y=lat, x=long+2*reef_offset, group=group), fill=clrs.cyclones[4],show.legend = FALSE) +
    geom_text(data=tops %>% mutate(min=min+2*reef_offset, max=max+2*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Ita (Apr 2014)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+2*reef_offset, group=group),fill=NA,color='black', size=0.2) +      
    geom_path(data=ita.fort %>% filter(long<150 | lat < -17), aes(y=lat, x=long+2*reef_offset), color='blue') +
    
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+3*reef_offset, group=group),fill='grey', color='grey70') +    
    geom_polygon(data=marcia.4m.fort, aes(y=lat, x=long+3*reef_offset, group=group), fill=clrs.cyclones[2],show.legend = FALSE) +
    geom_polygon(data=marcia.6m.fort, aes(y=lat, x=long+3*reef_offset, group=group), fill=clrs.cyclones[3],show.legend = FALSE) +
    geom_polygon(data=marcia.8m.fort, aes(y=lat, x=long+3*reef_offset, group=group), fill=clrs.cyclones[4],show.legend = FALSE) +
    geom_text(data=tops %>% mutate(min=min+3*reef_offset, max=max+3*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Marcia (Feb 2015)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+3*reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_path(data=marcia.fort %>% filter(long<153), aes(y=lat, x=long+3*reef_offset), color='blue') +

    
    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+4*reef_offset, group=group),fill='grey', color='grey70') +
    geom_polygon(data=nathan.4m.fort, aes(y=lat, x=long+4*reef_offset, group=group, fill='Hs 4')) +
    geom_polygon(data=nathan.6m.fort, aes(y=lat, x=long+4*reef_offset, group=group, fill='Hs 6')) +
    geom_polygon(data=nathan.8m.fort, aes(y=lat, x=long+4*reef_offset, group=group, fill='Hs 8')) +
    geom_text(data=tops %>% mutate(min=min+4*reef_offset, max=max+4*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Nathan (Mar 2015)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+4*reef_offset, group=group),fill=NA,color='black', size=0.2) +  
    geom_path(data=nathan.fort, aes(y=lat, x=long+4*reef_offset), color='blue') +

    geom_polygon(data=gbr.fort %>% filter(id=='Reef'), aes(y=lat, x=long+5*reef_offset, group=group),fill='grey', color='grey70') +
    geom_polygon(data=debbie.4m.fort, aes(y=lat, x=long+5*reef_offset, group=group, fill='Hs 4')) +
    geom_polygon(data=debbie.6m.fort, aes(y=lat, x=long+5*reef_offset, group=group, fill='Hs 6')) +
    geom_polygon(data=debbie.8m.fort, aes(y=lat, x=long+5*reef_offset, group=group, fill='Hs 8')) +
    geom_text(data=tops %>% mutate(min=min+5*reef_offset, max=max+5*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='Debbie (Mar 2017)'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+5*reef_offset, group=group),fill=NA,color='black', size=0.2) +  
    geom_path(data=debbie.fort, aes(y=lat, x=long+5*reef_offset), color='blue') +
    
    scale_fill_manual('', breaks = c('Hs 0','Hs 4', 'Hs 6', 'Hs 8'), labels=c(expression(Hs<4*m),expression(Hs>=4*m~'&'<6*m), expression(Hs>=6*m~'&'<8*m), expression(Hs>=8*m)), values=clrs.cyclones, limits=c('Hs 0','Hs 4', 'Hs 6', 'Hs 8'))+
    scale_x_continuous(breaks=seq(142,154,by=2), position = 'top') +
    scale_y_continuous(breaks=seq(-25,-10,by=2)) +
    coord_equal(ylim=bbox(management)[2,]+c(0,1),xlim=c(bbox(management)[1,]-1,bbox(management)[1,]+39)) +
    theme_minimal() +
    theme(legend.position=c(0.01,0.01), legend.justification=c(0,0),
          legend.background = element_blank(),
          legend.text.align = 0,legend.key=element_rect(color='black',fill=NA),
          axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          panel.border = element_rect(color='black',fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
#g2+    annotation_custom(grob=legd.cyclones,xmin=-Inf,ymin = -Inf,xmax=147,ymax=-21)


## Put the two figures together
#g1a <- ggplotGrob(g1+ theme(plot.margin=unit(c(0,0,0,0),'lines'), panel.spacing = unit(0,'lines'))
g1a <- ggplotGrob(g1+ theme(plot.margin=margin(t=0,b=-0.3,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')))
g2a <- ggplotGrob(g2+ theme(plot.margin=margin(t=0,b=0,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')))
g <- rbind(g1a, g2a, size="first") # stack the two plots
g$widths <- unit.pmax(g1a$widths, g2a$widths) # use the largest widths
grid.newpage()
grid.draw(g)

ggsave(filename=paste0('output/figures/Fig1mapC_',ZONES,'Zones_new.pdf'), g, width=10,height=7.3)
#ggsave(filename=paste0('output/figures/Fig1mapC_',ZONES,'Zones.pdf'), g, width=10,height=7.3)




## Lets see about adding a third panel for the coral change....
## I realize that something similar is calculated earlier, but I am now going to do the calculations again...
load('data/processed/manta.sum.RData')
manta.change=manta.sum %>% filter(REPORT_YEAR>=2010) %>%
    group_by(Location,Zone,REEF_NAME,Longitude,Latitude) %>%
    mutate(Cover.lag = lag(Cover),
           Change = 100*(Cover-Cover.lag),#/Cover.lag,
           Time.lag=lag(REPORT_YEAR),
           Time.change=REPORT_YEAR-Time.lag,
           Change.annual=Change/Time.change) %>%
    filter(REPORT_YEAR>2012) %>%
    ungroup

clrs.coral = c('#FF0000','#FF00FF')
legd.coral <- legendGrob(c("Decrease","Increase"), pch=21,gp=gpar(col = 'black', fill = clrs.coral),vgap=unit(0.15,'lines'))
#colLim=levels(coral.2016$Cat)
lgnd.dat = data.frame(x=c(149,149),y=c(-16,-16.7)) #ff(strt,coefs,length=2)
lgnd.dat$lab = c("Decrease","Increase")
lgnd.dat$Cat=factor(c(FALSE,TRUE))
lgnd.dat.size = ff(lgnd.dat[2,1:2],coefs,length=5)
lgnd.dat.size = ff(strt,coefs,length=5)
lgnd.dat.size$lab = format(seq(2.5,12.5,by=2.5),digits=2)
lgnd.dat.size$size = seq(2.5,12.5,by=2.5)

reef_offset = 8
g3 = ggplot(data=NULL, aes(y=Latitude,x=Longitude)) +
    geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2013) %>% mutate(Longitude=Longitude),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min, max=max, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2013'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2014) %>% mutate(Longitude=Longitude+reef_offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min+reef_offset, max=max+reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2014'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+2*reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2015) %>% mutate(Longitude=Longitude+2*reef_offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min+2*reef_offset, max=max+2*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2015'), vjust=0)  +
    geom_polygon(data=management.df, aes(y=lat, x=long+3*reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2016) %>% mutate(Longitude=Longitude+3*reef_offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min+3*reef_offset, max=max+3*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2016'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+4*reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2017) %>% mutate(Longitude=Longitude+4*reef_offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min+4*reef_offset, max=max+4*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2017'), vjust=0) +
    geom_polygon(data=management.df, aes(y=lat, x=long+5*reef_offset, group=group),fill=NA,color='black', size=0.2) +
    geom_point(data=manta.change %>% filter(REPORT_YEAR==2018) %>% mutate(Longitude=Longitude+5*reef_offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
    geom_text(data=tops %>% mutate(min=min+5*reef_offset, max=max+5*reef_offset, long=mean(c(min,max))), aes(y=lat+0.2, x=long,label='2018'), vjust=0) +
    scale_fill_manual('',breaks=c(FALSE,TRUE), values=c('red','green'), labels=c('Decrease','Increase')) +
    #scale_size_discrete('Annual change (%)') +
    #geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
    #geom_text(data=lgnd.dat, aes(y=y,x=x+0.3+0.3, label=lab), size=3,hjust=0, parse=TRUE) +
    #geom_point(data=lgnd.dat, aes(y=y,x=x+0.3, fill=Cat), shape=21, size=3) +
    #geom_point(data=lgnd.dat.size, aes(y=y,x=x+0.4, size=size), shape=21) +
    #geom_text(data=lgnd.dat.size, aes(y=y,x=x+0.4+0.4, label=lab), size=3,hjust=0, parse=TRUE) +
    scale_x_continuous(breaks=seq(142,154,by=2), position = 'top') +
    scale_y_continuous(breaks=seq(-25,-10,by=2)) +
    coord_equal(ylim=bbox(management)[2,]+c(0,1),xlim=c(bbox(management)[1,]-1,bbox(management)[1,]+39)) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          panel.border = element_rect(color='black',fill=NA),
          legend.position=c(0.99,0.99), legend.justification=c(0,1),
          #legend.position=c(0.9,0.5), legend.justification=c(1,1), legend.box='horizontal', legend.margin=margin(r=50,l=50),
          legend.background = element_blank(),
          legend.text.align = 0,
          plot.margin=margin(t=0,b=0,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')
          #legend.key=element_rect(color='black',fill=NA),
          )    +
    scale_radius('Annual\nchange (%)', breaks=c(5,10,15,20,25)) +
    guides(fill=guide_legend(order=1, label.theme=element_text(size=10), override.aes=list(size=5)),
           size=guide_legend(order=0))

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

legend <- g_legend(g3)

plotNew <- g3 + 
    #annotation_custom(grob = legend[2,], xmin = 143, xmax = 145, ymin = -19, ymax = -23) +
                                        #annotation_custom(grob = legend[-2,], xmin = 191, xmax = 195, ymin =tops$lat, ymax = -15.7)
    annotation_custom(grob = legend[-2,], xmin = 141, xmax = 145, ymin = -18, ymax = -21) +
    annotation_custom(grob = legend[2,], xmin = 190, xmax = 193, ymin =tops$lat, ymax = -15.7)
tmp <- ggplot_gtable(ggplot_build(plotNew))
tmp=cowplot::gtable_remove_grobs(tmp, "guide-box")
grid.newpage()
grid.draw(tmp)

g3a=tmp


## Put the three figures together
#g1a <- ggplotGrob(g1+ theme(plot.margin=unit(c(0,0,0,0),'lines'), panel.spacing = unit(0,'lines'))
g1a <- ggplotGrob(g1+ theme(plot.margin=margin(t=0,b=-0.3,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')))
g2a <- ggplotGrob(g2+ theme(plot.margin=margin(t=0,b=-0.3,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')))
#g3a <- ggplotGrob(g3+ theme(plot.margin=margin(t=0,b=0,l=0,r=0,unit='lines'), panel.spacing = unit(0,'lines')))
g <- rbind(g1a, g2a, g3a, size="first") # stack the three plots
g$widths <- unit.pmax(g1a$widths, g2a$widths, g3a$widths) # use the largest widths
grid.newpage()
grid.draw(g)

ggsave(filename=paste0('output/figures/Fig1mapE_',ZONES,'Zones_new.pdf'), g, width=10,height=10.2)



## Now for a figure that has the annual change for each year.
load('data/processed/manta.sum.RData')
manta.change=manta.sum  %>%
    group_by(Location,Zone,REEF_NAME,Longitude,Latitude) %>%
    mutate(Cover.lag = lag(Cover),
           Change = 100*(Cover-Cover.lag),#/Cover.lag,
           Time.lag=lag(REPORT_YEAR),
           Time.change=REPORT_YEAR-Time.lag,
           Change.annual=Change/Time.change) %>%
    filter(REPORT_YEAR>min(REPORT_YEAR)) %>%
    ungroup %>%
    droplevels

gA=ggplot(data=NULL, aes(y=Latitude,x=Longitude)) #+
    #geom_polygon(data=management.df, aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2)
long.offset=8
lat.offset=-17
Yrs = sort(unique(manta.change$REPORT_YEAR))
NROW=8
NCOL=length(Yrs)/NROW
mat=matrix(c(sort(unique(manta.change$REPORT_YEAR))), ncol=8, byrow=TRUE)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        gA=gA +
            geom_polygon(data=management.df %>% mutate(lat=lat+(i-1)*lat.offset, long=long+(j-1)*long.offset), aes(y=lat, x=long, group=group),fill=NA,color='black', size=0.2) +
            geom_point(data=manta.change %>% filter(REPORT_YEAR==mat[i,j]) %>% mutate(Latitude=Latitude+(i-1)*lat.offset, Longitude=Longitude+(j-1)*long.offset),aes(fill=Change.annual>0, size=abs(Change.annual)),shape=21,alpha=0.4) +
            geom_text(data=tops %>% mutate(lat=lat+0.2 + (i-1)*lat.offset, min=min+(j-1)*long.offset, max=max+(j-1)*long.offset, long=mean(c(min,max))), aes(y=lat, x=long),label=mat[i,j], vjust=0) 
    }
}

gA=gA + coord_equal() +
    scale_fill_manual('',breaks=c(FALSE,TRUE), values=c('red','green'), labels=c('Decrease','Increase')) +
    scale_radius('Annual\nchange (%)', breaks=c(5,10,15,20,25)) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    guides(fill=guide_legend(order=1, label.theme=element_text(size=10), override.aes=list(size=3)),
           size=guide_legend(order=0, label.theme=element_text(size=10)))

ggsave(filename=paste0('output/figures/AnnualChange_',ZONES,'Zones_new.pdf'), gA, width=10,height=10.2)
save(manta.change, file='data/processed/manta.change.RData')
