source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

load('data/processed/manta.sum.RData')
load(file='data/modelled/dat.gbr.RData')
load(file='data/modelled/last_year.RData')
load(file='data/modelled/mod.gbr.RData')




CoralTrends_changeInYears = function(mod,yr1,yr2) {
    load('data/processed/manta.sum.RData')    
    a=fixef(mod)[1]
    b=fixef(mod)[-1]
    a.reef=ranef(mod)[[1]][,'(Intercept)']
    b.reef=as.matrix(ranef(mod)[[1]][,-1])
    
    bb=sweep(b.reef,2,b,FUN='+')
    bbb = sweep(bb,1,(a+a.reef), FUN='+')
    state = data.frame(binomial()$linkinv((bbb)))
    state$REEF_NAME = rownames(state)
    state = state %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))
    
    last_year = state %>% dplyr:::select_('REEF_NAME',paste0('Year',yr2),paste0('Year',yr1)) %>%
        mutate_(.dots=setNames(paste0('Year',yr2,'-Year',yr1),'Diff')) %>%
        mutate(D=Diff<0) %>% dplyr:::select(-starts_with('Year'))
    last_year = last_year %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))
    last_year
}

CoralTrends_predict_reef_year = function(mod) {
    load('data/processed/manta.sum.RData')    
    a=fixef(mod)[1]
    b=fixef(mod)[-1]
    a.reef=ranef(mod)[[1]][,'(Intercept)']
    b.reef=as.matrix(ranef(mod)[[1]][,-1])
    
    bb=sweep(b.reef,2,b,FUN='+')
    bbb = sweep(bb,1,(a+a.reef), FUN='+')
    state = data.frame(binomial()$linkinv((bbb)))
    state$REEF_NAME = rownames(state)
    state = state %>% left_join(manta.sum %>% group_by(REEF_NAME) %>% summarize_at(vars(Latitude,Longitude,Tows), funs(mean)))
    state
}

load(file='data/modelled/mod.northern.RData')
last_year.northern = CoralTrends_changeInYears(mod.northern,yr1='2015',yr2='2017')

load(file='data/modelled/mod.central.RData')
last_year.central = CoralTrends_changeInYears(mod.central,yr1='2016',yr2='2017')

load(file='data/modelled/mod.southern.RData')
last_year.southern = CoralTrends_changeInYears(mod.southern,yr1='2016',yr2='2017')

last_year = rbind(last_year.northern, last_year.central, last_year.southern)

library(ggsn)

ewbrks <- seq(144,152,by=2)
nsbrks <- seq(-24,-10,by=2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste0(x, "°W"), ifelse(x > 0, paste0(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), "°S"), ifelse(x > 0, paste0(x, "°N"),x))))

gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year %>% arrange(desc(D)), aes(y=Latitude,x=Longitude,fill=D), size=2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, dd2km = TRUE, model = 'WGS84',st.size=3,location='bottomleft') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
ggsave(file='output/figures/InfluenceMap_manta.png', gp, width=5, height=5, dpi=300)
ggsave(file='output/figures/InfluenceMap_manta.pdf', gp, width=5, height=5, dpi=300)


gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=abs(last_year$Diff)*last_year$Tows,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, dd2km = TRUE, model = 'WGS84',st.size=3,location='bottomleft') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
gp
ggsave(file='output/figures/InfluenceMap_mantaBubble.png', gp, width=5, height=5, dpi=300)
ggsave(file='output/figures/InfluenceMap_mantaBubble.pdf', gp, width=5, height=5, dpi=300)


## Now we will repeat this, yet only plot the reefs that were actually sampled in the last year
manta.sum.reefs = manta.sum %>% filter(REPORT_YEAR==2017) %>% dplyr:::select(REEF_NAME) %>% distinct
gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=last_year %>% right_join(manta.sum.reefs), aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,size=abs(last_year$Diff)*last_year$Tows), shape=21, alpha=0, color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    #geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=last_year %>% right_join(manta.sum.reefs) %>% arrange(desc(D)), aes(y=Latitude,x=Longitude,fill=D), size=2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    #scale_size_area('Influence', max_size=4)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.75),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, dd2km = TRUE, model = 'WGS84',st.size=3,location='bottomleft') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
ggsave(file='output/figures/InfluenceMap_manta_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='output/figures/InfluenceMap_manta_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)


ly=last_year %>% right_join(manta.sum.reefs)
gp=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=Diff<0), size=0,shape=21, alpha=0, color='black') +
    geom_point(data=ly, aes(y=Latitude,x=Longitude,size=abs(Diff)*100), shape=21, alpha=0, color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0), size=abs(100*last_year$Diff), alpha=0.5,shape=21,color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=Diff<0),  size=last_year$Tows/20,alpha=0.5,shape=21,color='black') +
    ##geom_point(data=last_year, aes(y=Latitude,x=Longitude,fill=D), size=(abs(last_year$Diff)*last_year$Tows)/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    #geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D), size=abs(ly$Diff)*ly$Tows,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    geom_point(data=ly, aes(y=Latitude,x=Longitude,fill=D), size=abs(ly$Diff)*100/2,alpha=0.5,shape=21,color='black', show.legend = FALSE) +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('', breaks = ewbrks, labels = ewlbls) +#, expand = c(0, 0)) +
    scale_y_continuous('',breaks = nsbrks, labels = nslbls) + #, expand = c(0, 0))+
    scale_fill_manual('Change', breaks=c(FALSE,TRUE), labels=c('Increase','Decrease'), values=c('blue','red'),limits=c(FALSE,TRUE))+
    ##scale_size('Magnitude of change', breaks=c(4,8,12,16), labels=c(4,8,12,16),values=c(4,8,12,16))+
                                        #scale_size_identity('Magnitude of change') +
    scale_size('Magnitude', breaks=c(0.5,1,2,5,10,15),range=c(0.5,15)/2) +
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            #legend.position = c(0.95,0.95),legend.justification=c(1,1),
                            legend.position = c(0.95,0.65),legend.justification=c(1,0.5),legend.background = element_rect(fill='white', color=NA))+
    guides(fill=guide_legend(override.aes = list(size=5,alpha=0.5)),
           size=guide_legend(override.aes = list(alpha=0.5)),
           color=guide_legend(override.aes = list(size=5,alpha=1)))
gp=gp + scalebar(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2],
                 dist=200, dd2km = TRUE, model = 'WGS84',st.size=3,location='bottomleft') +
    north(x.min=bbox(whagbr)[1,1], x.max=bbox(whagbr)[1,2], y.min=bbox(whagbr)[2,1],y.max=bbox(whagbr)[2,2], location='topright',scale=0.1,symbol=12)
gp

ggsave(file='output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.png', gp, width=7, height=7, dpi=300)
ggsave(file='output/figures/InfluenceMap_mantaBubble_OnlyLastYearReefs.pdf', gp, width=7, height=7, dpi=300)



    
