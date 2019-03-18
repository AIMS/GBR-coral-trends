source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

INCLUDE_GBR = FALSE

load('data/processed/manta.sum.RData')
load(file='data/modelled/dat.gbr.RData')
load(file='data/modelled/dat.northern.RData')
load(file='data/modelled/dat.central.RData')
load(file='data/modelled/dat.southern.RData')

load(file='data/spatial/spatial_3Zone.RData')

load(file='data/modelled/cots.sum.all_3Zone.RData')
load(file='data/modelled/bleaching.sum.all_3Zone.RData')
load(file='data/modelled/cyclones.sum.all_3Zone.RData')


## We have decided to bring the GBR wide figure out on its own
newdata=dat.gbr %>% mutate(Location=factor(Location,labels=c('Great Barrier\n\nReef')))
max_year = max(as.numeric(as.character(newdata$Year)))
nd = newdata %>%
    group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")
g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    annotate(geom='rect', xmin=2012, xmax=max_year, ymin=-Inf,ymax=Inf, fill='grey70',alpha=0.4) +
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,35)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,35)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2020))+
                theme_classic()+
                    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                          panel.background=element_rect(color='black'),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color=NA),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                          strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          plot.margin=unit(c(0,0,2,0),'pt'))
## Now we generate the banner
a=oz:::ozRegion(sections=c(3,11:13))
a=oz:::ozRegion()
cc=rbind(xy2df(a$lines[[3]]),
    xy2df(a$lines[[13]]),
    xy2df(a$lines[[12]])[nrow(xy2df(a$lines[[12]])):1,],
    xy2df(a$lines[[11]]))
aa.ps<-SpatialPolygons(list(Polygons(list(Polygon(cc)),ID="QLD")))

gt = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    #geom_blank(aes(x=270,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  
gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
gt0=gt1

gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg.gbr <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=4, b=b, r=4, name="pic_predator"))

grid.draw(gg.gbr)
save(gg.gbr, file='data/spatial/gg.gbr_3Zone.RData')
ggsave(file='output/figures/GBRStan_3Zone.pdf', grid.draw(gg.gbr), width=4, height=3, units='in',dpi=300) 
ggsave(file='output/figures/GBRStan_3Zone.png', grid.draw(gg.gbr), width=4, height=3, units='in',dpi=300) 



## Just whole GBR
newdata = rbind(dat.gbr,dat.northern,dat.central,dat.southern)
newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location),
                           labels=c('Great Barrier\n\nReef','Northern\n', 'Central\n', 'Southern\n'))
newdata = newdata %>% filter(Location=='Great Barrier\n\nReef') %>% droplevels
nd = newdata %>%
    group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")
g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2020))+
                theme_classic()+
                    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                          panel.background=element_rect(color='black'),
                          #axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=0.5,unit='lines')),
                          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color=NA),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                          #strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          plot.margin=unit(c(0,0,2,0),'pt'))
a=oz:::ozRegion(sections=c(3,11:13))
a=oz:::ozRegion()
cc=rbind(xy2df(a$lines[[3]]),
    xy2df(a$lines[[13]]),
    xy2df(a$lines[[12]])[nrow(xy2df(a$lines[[12]])):1,],
    xy2df(a$lines[[11]]))
aa.ps<-SpatialPolygons(list(Polygons(list(Polygon(cc)),ID="QLD")))

gt = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    #geom_blank(aes(x=270,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  

gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4]) #+
    #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)
## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
grid.draw(gg)
### The following (and above) produces the three panels with the fancy facets 
save(gg, file='data/spatial/gg_3Zone_GBR.RData')
ggsave(file='output/figures/3ZonesGBR.pdf', grid.draw(gg), width=5, height=5, units='in',dpi=300) 
ggsave(file='output/figures/3ZonesGBR.png', grid.draw(gg), width=5, height=5, units='in',dpi=300)




#Start with all panel plot============================================================================================
newdata = rbind(dat.gbr,dat.northern,dat.central,dat.southern)
newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location),
                           labels=c('Great Barrier\n\nReef','Northern\n', 'Central\n', 'Southern\n'))
save(newdata, file='data/modelled/cellmeans.RData')
if (!INCLUDE_GBR) newdata = newdata %>% filter(Location!='Great Barrier\n\nReef') %>% droplevels
#newdata = newdata %>% left_join(dat.all %>% select(Location, N) %>% distinct)
nd = newdata %>%
    group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")
g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2020))+
                theme_classic()+
                    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                          panel.background=element_rect(color='black'),
                          #axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=0.5,unit='lines')),
                          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color=NA),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                          #strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          plot.margin=unit(c(0,0,2,0),'pt'))

## Get the dimensions of the strips
#gg=ggplotGrob(g1)
#index <- which(sapply(gg$grobs, function(x) x$name == "strip"))
#gg$grobs[[index[1]]]$heights
#gg$grobs[[index[1]]]$widths = gg$grobs[[index[1]]]$heights * 2
#grid.draw(gg)
#grep("strip-t-1-1", gg$layout$name)

## Now we generate the banner
a=oz:::ozRegion(sections=c(3,11:13))
a=oz:::ozRegion()
cc=rbind(xy2df(a$lines[[3]]),
    xy2df(a$lines[[13]]),
    xy2df(a$lines[[12]])[nrow(xy2df(a$lines[[12]])):1,],
    xy2df(a$lines[[11]]))
aa.ps<-SpatialPolygons(list(Polygons(list(Polygon(cc)),ID="QLD")))

gt = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    #geom_blank(aes(x=270,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  

gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4]) #+
    #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)

gt2=gt+geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt3=gt+geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])

gt4=gt+geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])


if (!INCLUDE_GBR) gt1=gt2; gt2=gt3; gt3=gt4; gt4=gt5;

## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
#facets <- grep("strip-t-4-1", gT$layout$name)
#gg <- with(gg$layout[facets,],
#           gtable_add_grob(gg, ggplotGrob(gt4),t=t, l=16, b=b, r=16, name="pic_predator"))
if (INCLUDE_GBR) {
    facets <- grep("strip-t-5-1", gT$layout$name)
                                        #gg <- with(gg$layout[facets,],
                                        #           gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=20, b=b, r=20, name="pic_predator"))
    gg <- with(gg$layout[facets,],
               gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=17, b=b, r=20, name="pic_predator"))
}

grid.draw(gg)
### The following (and above) produces the three panels with the fancy facets 
save(gg, file='data/spatial/gg_3Zone.RData')
#ggsave(file='output/figures/3ZonesStan.pdf', grid.draw(gg), width=15, height=3, units='in',dpi=300) 
#ggsave(file='output/figures/3ZonesStan.png', grid.draw(gg), width=15, height=3, units='in',dpi=300) 

ggsave(file='output/figures/3Zones.pdf', grid.draw(gg), width=10, height=3, units='in',dpi=300) 
ggsave(file='output/figures/3Zones.png', grid.draw(gg), width=10, height=3, units='in',dpi=300) 

#============================================================================================================

##gg$grobs[[31]]$heights
#gg$grobs[[31]]$widths
#gh = ggplotGrob(gt5)
#gh$grobs[[grep("panel", gh$layout$name)]]$widths


## Now the cots data
if (!INCLUDE_GBR) cots.sum.all = cots.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
labs = levels(cots.sum.all$Location)
#Version with NO
gcots = ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('NO','AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
    geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill='blue',show.legend=FALSE) +
    #geom_line(stat='identity',position='stack',aes(alpha=COTScat), fill='red',show.legend=FALSE) +
    #facet_grid(~Location) +
    #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(90,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    theme_classic() +
    theme(
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,0,0,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
#Version without NO
gcots = ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
    geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill='blue',show.legend=FALSE) +
    #geom_line(stat='identity',position='stack',aes(alpha=COTScat), fill='red',show.legend=FALSE) +
    #facet_grid(~Location) +
    #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    scale_alpha_manual(breaks=c('IO','AO'), values=c(0.4,1)) +
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(90,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    theme_classic() +
    theme(
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,0,0,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
gcots

gcots <- ggplot_gtable(ggplot_build(gcots))
gcots$widths = gg$widths

grid.arrange(gg,gcots, nrow=2,heights=c(2,1),padding=unit(-1,'lines'))



## Combine cots and bleaching
#csa = cots.sum.all %>% dplyr:::select(Zone,REPORT_YEAR,Location,COTScat,COTS.p) %>% gather(key='Impact', value='Percent', COTS.p)
#bsa = bleaching.sum.all %>% dplyr:::select(Zone,REPORT_YEAR,Location,BLEACHINGcat,BLEACHING.p) %>% gather(key='Impact', value='Percent', BLEACHING.p)

#gImpact = cots.sum.all %>% full_join(bleaching.sum.all)
library(viridis)
library(RColorBrewer)
labs = levels(cots.sum.all$Location)
cots.dat = cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('IO','AO'))
gcots = ggplot(cots.dat, aes(y=COTS.p, x=REPORT_YEAR-0.3)) +
    geom_bar(stat='identity',position='stack',aes(fill=COTScat),width=0.1,show.legend=FALSE)+
    geom_bar(data=cots.dat %>% filter(COTScat %in% c('IO','AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.4,show.legend=FALSE)+
    geom_bar(data=cots.dat %>% filter(COTScat %in% c('AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.4,show.legend=FALSE)+
    #geom_col(position='stack',aes(fill=COTScat),width=0.1,show.legend=TRUE)+
    #geom_area(position='stack',aes(fill=COTScat),show.legend=TRUE)+
    scale_fill_manual('', breaks=c('IO','AO'), labels=c('IO','AO'), values=scales:::brewer_pal(palette='Greens')(3)[-1]) +
    #geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill=brewer.pal(3,'Dark2')[1],show.legend=FALSE,width=0.3) +
    #scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    theme_classic() +
    theme(
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,0,5,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
gcots

if (!INCLUDE_GBR) bleaching.sum.all = bleaching.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
bleaching.dat = bleaching.sum.all %>% filter(REPORT_YEAR >1985, BLEACHINGcat!='0')
gbleaching = ggplot(bleaching.dat, aes(y=BLEACHING.p, x=REPORT_YEAR+0.3)) +
    geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('2','3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
    
    #geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),show.legend=FALSE,width=0.3)+
    scale_fill_manual('', breaks=c(1,2,3,4), labels=c(1,2,3,4), values=c(scales:::brewer_pal(palette='Reds')(5)[-1])) +
    #geom_bar(stat='identity',position='stack',aes(alpha=BLEACHINGcat), fill=brewer.pal(3,'Dark2')[2],show.legend=FALSE,width=0.3) +
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    #scale_alpha_manual(breaks=c(0,1,2,3,4), values=c(0,0.25,0.5,0.75,1)) +
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    theme_classic() +
    theme(
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,0,5,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
gbleaching

if (!INCLUDE_GBR) cyclones.sum.all = cyclones.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
cyclones.dat = cyclones.sum.all %>% filter(REPORT_YEAR >1985, CYCLONEcat!=0)
gcyclones = ggplot(cyclones.dat, aes(y=CYCLONE.p, x=REPORT_YEAR+0)) +
    geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+
    geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('2','3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+
    geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+

                                        #geom_bar(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill=brewer.pal(3,'Dark2')[3],show.legend=FALSE,width=0.3) +
    #geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),show.legend=FALSE,width=0.3)+
    scale_fill_manual('', breaks=c(1,2,3), labels=c(1,2,3), values=c(scales:::brewer_pal(palette='Blues')(4)[-1])) +
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    #scale_alpha_manual(breaks=c(0,1,2,3,4,5), values=c(0,1/5,2/5,3/5,4/5,5/5)) +
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    theme_classic() +
    theme(
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,0,5,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
gcyclones


gbleaching = gbleaching + theme(panel.background=element_blank())
gcyclones = gcyclones + theme(panel.background=element_blank())
gcots <- ggplot_gtable(ggplot_build(gcots))
gbleaching <- ggplot_gtable(ggplot_build(gbleaching))
gcyclones <- ggplot_gtable(ggplot_build(gcyclones))
panels <- grepl("panel", gbleaching$layout$name)

pp <- c(subset(gcots$layout, grepl("panel", gcots$layout$name), se = t:r))

# Overlap panels for second plot on those of the first plot
gT <- gtable_add_grob(gcots, gbleaching$grobs[grepl("panel", gcots$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l, name='bleaching')
gT <- gtable_add_grob(gT, gcyclones$grobs[grepl("panel", gcots$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l, name='cyclones')
gT$widths = gg$widths

grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines'))
save(gT, file='data/spatial/gT_3Zone.RData')
ggsave(file='output/figures/3ZonesFigurePts3_new.pdf', grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines')), width=10, height=5, units='in',dpi=300) 
ggsave(file='output/figures/3ZonesFigurePts1_new.png', grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines')), width=15, height=5, units='in',dpi=300) 

## Stop here - not quite - just do the next bit..



## As an alternative, make a stand alone version of the disturbances that can be added to the modelled disturbances
## NOTE, TO ADJUST THE WIDTH AND POSITION OF THE BARS, IT IS NECESSARY TO ALTER THE x POSITION AND WIDTH 
library(viridis)
library(RColorBrewer)
if (!INCLUDE_GBR) cots.sum.all = cots.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
labs = levels(cots.sum.all$Location)
cots.dat = cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('IO','AO'))
gcots = ggplot(cots.dat, aes(y=COTS.p, x=REPORT_YEAR-0.3)) +
    geom_bar(stat='identity',position='stack',aes(fill=COTScat),width=0.1,show.legend=FALSE)+
    geom_bar(data=cots.dat %>% filter(COTScat %in% c('IO','AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.3,show.legend=FALSE)+
    geom_bar(data=cots.dat %>% filter(COTScat %in% c('AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.3,show.legend=FALSE)+
    #geom_col(position='stack',aes(fill=COTScat),width=0.1,show.legend=TRUE)+
    #geom_area(position='stack',aes(fill=COTScat),show.legend=TRUE)+
    scale_fill_manual('', breaks=c('IO','AO'), labels=c('IO','AO'), values=scales:::brewer_pal(palette='Greens')(3)[-1]) +
    #geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill=brewer.pal(3,'Dark2')[1],show.legend=FALSE,width=0.3) +
    #scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
    facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    theme_classic() +
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,5,5,0),'pt'),
        strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
        #strip.text.x=element_blank(), strip.background=element_blank())
gcots

if (!INCLUDE_GBR) bleaching.sum.all = bleaching.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
bleaching.dat = bleaching.sum.all %>% filter(REPORT_YEAR >1985, BLEACHINGcat!='0')
gbleaching = ggplot(bleaching.dat, aes(y=BLEACHING.p, x=REPORT_YEAR+0.3)) +
    geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('2','3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
    geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
    #geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),show.legend=FALSE,width=0.3)+
    scale_fill_manual('', breaks=c(1,2,3,4), labels=c(1,2,3,4), values=c(scales:::brewer_pal(palette='Reds')(5)[-1])) +
    #geom_bar(stat='identity',position='stack',aes(alpha=BLEACHINGcat), fill=brewer.pal(3,'Dark2')[2],show.legend=FALSE,width=0.3) +
    facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    #scale_alpha_manual(breaks=c(0,1,2,3,4), values=c(0,0.25,0.5,0.75,1)) +
                                        #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    theme_classic() +
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,5,5,0),'pt'),
        strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
        #strip.text.x=element_blank(), strip.background=element_blank())
gbleaching

if (!INCLUDE_GBR) cyclones.sum.all = cyclones.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
cyclones.dat = cyclones.sum.all %>% filter(REPORT_YEAR >1985, CYCLONEcat!=0)
gcyclones = ggplot(cyclones.dat, aes(y=CYCLONE.p, x=REPORT_YEAR+0)) +
    geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=FALSE)+
    geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('2','3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=FALSE)+
    geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=FALSE)+
                                        #geom_bar(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill=brewer.pal(3,'Dark2')[3],show.legend=FALSE,width=0.3) +
    #geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),show.legend=FALSE,width=0.3)+
    scale_fill_manual('', breaks=c(1,2,3), labels=c(1,2,3), values=c(scales:::brewer_pal(palette='Blues')(4)[-1])) +
    facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    #scale_alpha_manual(breaks=c(0,1,2,3,4,5), values=c(0,1/5,2/5,3/5,4/5,5/5)) +

    scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
    #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    theme_classic() +
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,5,5,0),'pt'),
        ##strip.text.x=element_blank(), strip.background=element_blank())
        strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
gcyclones

gbleaching = gbleaching + theme(panel.background=element_blank()) + ggtitle('a)')
gcyclones = gcyclones + theme(panel.background=element_blank()) + ggtitle('a)')
gcots <- ggplot_gtable(ggplot_build(gcots + ggtitle('a)')))
gbleaching <- ggplot_gtable(ggplot_build(gbleaching))
gcyclones <- ggplot_gtable(ggplot_build(gcyclones))
panels <- grepl("panel", gbleaching$layout$name)

pp <- c(subset(gcots$layout, grepl("panel", gcots$layout$name), se = t:r))

# Overlap panels for second plot on those of the first plot
gT <- gtable_add_grob(gcots, gbleaching$grobs[grepl("panel", gcots$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l, name='bleaching')
gT <- gtable_add_grob(gT, gcyclones$grobs[grepl("panel", gcots$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l, name='cyclones')
#gT$widths = gg$widths

grid.draw(gT)

facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
grid.draw(gg)
save(gg, file='data/processed/disturbanceBars-gg.RData')
#grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines'))





## Annual change =================================================================================================
load('data/processed/manta.sum.RData')
## raw data
manta.change=manta.sum %>%
    group_by(REEF_NAME) %>%
    mutate(Cover.lag = lag(Cover),
           Change = 100*(Cover-Cover.lag),#/Cover.lag,
           Time.lag=lag(REPORT_YEAR),
           Time.change=REPORT_YEAR-Time.lag,
           Change.annual=Change/Time.change)

manta.change = manta.change %>% ungroup %>%
    arrange(REEF_ID,REPORT_YEAR) %>%
    mutate(Location=factor(Location, levels=unique(Location), labels=c('Northern\n', 'Central\n', 'Southern\n')))

sampling_years = manta.change %>% dplyr::select(Location,Zone,REEF_NAME,REPORT_YEAR) %>% distinct %>%
    arrange(Location,Zone,REEF_NAME,REPORT_YEAR)

## Fill in the years
manta.change =
    manta.change %>%
    group_by(REEF_NAME) %>%
    do({
        x=.
        x=x%>% full_join(
                   data.frame(Location=unique(x$Location),
                              Zone=unique(x$Zone),
                              REEF_NAME=unique(x$REEF_NAME),
                              REPORT_YEAR=seq(min(x$REPORT_YEAR),max(x$REPORT_YEAR), by=1))) %>%
            arrange(REPORT_YEAR)
        value = x$Change.annual
        for (i in nrow(x):1) {
            value[i] = ifelse(is.na(value[i]), value[i+1], value[i])
        }
        x = x %>% mutate(Change.annual=value)
        x
    }) %>% ungroup %>%
    mutate(Location=factor(Location, levels=unique(Location), labels=c('Northern\n', 'Central\n', 'Southern\n')))

## Remove the non-sampling years
manta.change =
    manta.change %>% mutate(Exclude=ifelse(Zone=='Northern' & REPORT_YEAR %in% c(2012,2014,2016), 'Y',NA)) %>%
    filter(is.na(Exclude)) %>%
    dplyr::select(-Exclude)
save(manta.change, file='data/modelled/manta.change.RData')

g2 <- ggplot(manta.change, aes(y=Change.annual, x=as.numeric(as.character(REPORT_YEAR)))) +
    geom_hline(yintercept=0, color='red') +
    geom_point(size=0.5) +
    facet_grid(~Location) +
    scale_y_continuous('Change in \n% coral cover') +
    scale_x_continuous('',breaks=seq(1985,2015,by=5), limits=c(1985,2020), expand=c(0,0),position='top') +
    theme_classic()+
    theme(strip.background=element_blank(), strip.text.x=element_blank(),
          panel.border=element_rect(fill=NA,color='black'),
                  axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,10,0,0),'pt'))

old=manta.change

manta.change.bars = manta.change %>% group_by(Location, REPORT_YEAR) %>%
    summarize(N=n(),
              N.p=100*sum(Change.annual>0 & Change.annual<=5)/N,
              N.p5=100*sum(Change.annual>5 & Change.annual<=10)/N,
              N.n=-100*sum(Change.annual<0 & Change.annual>= -5)/N,
              N.n5=-100*sum(Change.annual< -5 & Change.annual>= -10)/N,
              #N.p5=100*sum(Change.annual>10)/N,
              #N.n5=-100*sum(Change.annual< -5)/N,
              N.p10=100*sum(Change.annual>10)/N,
              N.n10=-100*sum(Change.annual< -10)/N
              ) %>%
    gather(key=Change, value=Percent,-Location,-REPORT_YEAR,-N) %>%
    #mutate(Change=factor(Change, levels=c('N.n10', 'N.n', 'N.p', 'N.p10')))
    mutate(Change=factor(Change, levels=c('N.n10', 'N.p10','N.n5', 'N.p5','N.n', 'N.p'))) 
save(manta.change.bars, file='data/modelled/manta.change.bars.RData')

#ggplot(manta.change.bars, aes(y=Percent, x=REPORT_YEAR)) + geom_bar(stat='Identity', aes(fill=Change), alpha=0.3) + facet_wrap(~Location) +
#    scale_fill_manual('Change', breaks=c('N.p10','N.p','N.n','N.n10'), values=c('darkred','darkgreen','red','green'))


load(file='data/modelled/dat.northern.diff.RData')
load(file='data/modelled/dat.central.diff.RData')
load(file='data/modelled/dat.southern.diff.RData')
annual = rbind(dat.northern.diff, dat.central.diff, dat.southern.diff)
annual$Location <- factor(annual$Location, levels=unique(annual$Location),
                          labels=c('Northern\n', 'Central\n', 'Southern\n'))
## Make the changes annual
annual = annual %>% mutate(Year0=stringr::str_sub(Years, start=8), Diff=as.integer(Year)-as.integer(Year0)) %>%
    mutate_at(vars(mean,lower,upper), function(x,Diff=.$Diff) x/Diff)
save(annual, file='data/modelled/annual.change.RData')

hues <- RColorBrewer::brewer.pal(4, "Blues")

g1 <-ggplot(data=manta.change.bars,aes(y=Percent, x=REPORT_YEAR, fill=Change)) +
    geom_bar(stat='Identity', alpha=1) +
    scale_fill_manual('Change', breaks=c('N.p10','N.p5','N.p','N.n','N.n5','N.n10'),
                      values=c('darkred','darkgreen','red','green','orange','lightgreen'),
                      labels=c('Increase > 10%','Increase 5-10%','Increase <5%','Decline < 5%', 'Decline 5-10%', 'Decline > 10%')) +
    facet_wrap(~Location, nrow=1, scales='fixed') +
        scale_y_continuous('Proportion of reefs') +
    scale_x_continuous('',breaks=seq(1985,2015,by=5), limits=c(1985,2020), expand=c(0,0)) +
        theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
          strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
          plot.margin=unit(c(0,0,2,0),'pt'))
## Add the little maps in the facet strips
## finally we put them together
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
facets <- grep("strip-t-2-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
facets <- grep("strip-t-3-1", gT$layout$name)
gg <- with(gg$layout[facets,],
           gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
#facets <- grep("strip-t-4-1", gT$layout$name)
#gg <- with(gg$layout[facets,],
#           gtable_add_grob(gg, ggplotGrob(gt4),t=t, l=16, b=b, r=16, name="pic_predator"))
if (INCLUDE_GBR) {
    facets <- grep("strip-t-5-1", gT$layout$name)
                                        #gg <- with(gg$layout[facets,],
                                        #           gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=20, b=b, r=20, name="pic_predator"))
    gg <- with(gg$layout[facets,],
               gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=17, b=b, r=20, name="pic_predator"))
}

## Blend the two figures better
p2 <- ggplot_gtable(ggplot_build(g2))
p1 <- gg
                                        #p2$widths[c(1,2,3,4)] = p1$widths[c(1,2,3,4)]
p1$widths[c(1,2,3,4)]=p2$widths[c(1,2,3,4)]
p2$widths[13]=p1$widths[17]+p1$widths[16]
grid.arrange(p1,p2, nrow=2,heights=c(2,1),padding=unit(-1,'lines'))
ggsave(file='output/figures/CoralChangeFig.pdf', grid.arrange(p1,p2, nrow=2,heights=c(2,1),padding=unit(-1,'lines')), width=10, height=5, units='in',dpi=300) 






