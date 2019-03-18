source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

load(file='data/primary/manta.RData')

pcode.mod = read.csv('data/primary/P_CODE_MOD.csv', strip.white=TRUE)
manta = manta %>% left_join(pcode.mod)

## Data are collected per Manta Tow.
## The most appropriate unit for these analyses is the reef level.
## - Only include data collected after 1985
## - Convert data into percent cover
## - Summarize the data to the reef/year level
manta.sum = manta %>%
    filter(REPORT_YEAR>1985) %>%
    mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%  
    group_by(P_CODE.mod,A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR) %>%
        summarise(Cover=mean(Cover, na.rm=TRUE), Tows=length(unique(TOW_SEQ_NO)), Latitude=mean(REEF_LAT, na.rm=TRUE), Longitude=mean(REEF_LONG, na.rm=TRUE)) %>%
            ungroup()


## Include only those reefs that have had more than 4 visits
wch <- table(manta.sum$REEF_NAME)
wch<-names(wch[wch>4])
manta.sum = manta.sum %>% filter(REEF_NAME %in% wch)


## Put the reefs into Locations (Zones)

#manta.sum = manta.sum %>% mutate(Location=CoralTrends_calc3ZoneLocations(Latitude))
manta.sum = CoralTrends_calc3ZoneLocations(manta.sum)
## which reefs are excluded...
#as.character(unique(filter(manta.sum, Location=='Outside')$REEF_NAME))
## Reefs exclude (Outside)
exc=manta.sum %>% filter(Location=='Outside')
manta.sum = manta.sum %>% filter(Location!='Outside') %>% droplevels %>% mutate(Zone=factor(Location, levels=c('Northern','Central','Southern')))
## number of reefs used
length(unique(manta.sum$REEF_NAME))
## number of reef surveys
nrow(manta.sum)
writeLines(paste0("NumberOfReefs=",length(unique(manta.sum$REEF_NAME)),
                  "\nNumberOfSurveys=",nrow(manta.sum)),
           con='data/processed/Manta.properties')
save(manta.sum, file='data/processed/manta.sum.RData')



##For comparison, we can also attempt to replicate what Glenn had done...
manta.sum.G = manta %>%
    filter(REPORT_YEAR>1985,REPORT_YEAR<2013) %>%
    mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%  
    group_by(A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR) %>%
        summarise(Cover=mean(Cover, na.rm=TRUE), Tows=length(unique(TOW_SEQ_NO)), Latitude=mean(REEF_LAT, na.rm=TRUE), Longitude=mean(REEF_LONG, na.rm=TRUE)) %>%
            ungroup()

wch <- table(manta.sum.G$REEF_NAME)
wch<-names(wch[wch>4])
#length(wch)  # number of reefs
manta.sum.G = manta.sum.G %>% filter(REEF_NAME %in% wch)
#dim(manta.sum) # number of reef surveys
manta.sum.G = manta.sum.G %>%
    CoralTrends_calc3ZoneLocations %>%
    mutate(Location=ifelse(Latitude>-11.8, 'Outside',as.character(Location)))
## Reefs exclude (Outside)
exc.G=manta.sum.G %>% filter(Location=='Outside')

manta.sum.G = manta.sum.G %>% filter(Location!='Outside')
nrow(manta.sum.G)
length(unique(manta.sum.G$REEF_NAME))

