source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

load('data/processed/manta.sum.RData')

## Generate a list of reefs we are using to help accumulate other sources of associated data
all.reefs = manta.sum %>%
    dplyr:::select(P_CODE.mod,REEF_NAME,REEF_ID,Latitude,Longitude) %>%
    group_by(REEF_NAME,REEF_ID) %>%
    summarize_at(vars(Latitude,Longitude), funs(mean)) %>%
    as.data.frame
write.csv(all.reefs,file='data/all.reefs_3Zone.csv', quote=FALSE, row.names=FALSE)
save(all.reefs, file='data/all.reefs.RData')

## Genuine stan cannot handle proportional data for binomial families (particularly when weights are applied)
## A work-around is to multiple the proportion by the weights and convert this into an integer 
dat.all = manta.sum %>% dplyr:::select(Cover, REEF_NAME, Tows,P_CODE.mod,Location,REPORT_YEAR) %>%
    mutate(Year=factor(REPORT_YEAR), N=length(unique(REEF_NAME))) %>% ungroup() %>%
    mutate(Cvr1 = as.integer(as.vector(Cover) * Tows), Cvr0 = Tows - Cvr1)


## GBR
dat.all.gbr = dat.all %>% droplevels
mod.gbr <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|P_CODE.mod/REEF_NAME), data=dat.all.gbr, family=binomial,iter=5000,warmup=2500,chains=3,cores=3)

dat.gbr = data.frame(Location='Great Barrier Reef',Year=unique(dat.all.gbr$Year), N=length(unique(dat.all.gbr$REEF_NAME))) 
Xmat = model.matrix(~Year, dat.gbr)

coefs = data.frame(mod.gbr) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
dat.gbr = cbind(dat.gbr,
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
)
save(dat.gbr, file='data/modelled/dat.gbr.RData')
save(dat.all.gbr, file='data/modelled/dat.all.gbr.RData')
save(mod.gbr, file='data/modelled/mod.gbr.RData')

l=levels(dat.all.gbr$Year)
last_year=ranef(mod.gbr)[[1]][,c('(Intercept)',paste0('Year',l[(length(l)-1):length(l)]))]
                                        #last_year=ranef(mod.gbr)[[1]][,c('Year2016','Year2017')]
ly = cbind(binomial()$linkinv(last_year[,1]-last_year[,2]),
                  binomial()$linkinv(last_year[,1]-last_year[,3])
                  )
last_year$Diff = ly[,1]-ly[,2]
last_year$REEF_NAME = rownames(last_year)
save(last_year, file='data/modelled/last_year.RData')


## Northern
dat.all.northern = dat.all %>% filter(Location=='Northern') %>% droplevels
#dat.all.northern$Cvr1 = as.integer(as.vector(dat.all.northern$Cover) * dat.all.northern$Tows)
#dat.all.northern$Cvr0 = dat.all.northern$Tows - dat.all.northern$Cvr1
mod.northern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|REEF_NAME), data=dat.all.northern, family=binomial,iter=5000,warmup=2500,chains=3,cores=3)

dat.northern = data.frame(Location='Northern',Year=unique(dat.all.northern$Year), N=length(unique(dat.all.northern$REEF_NAME))) 
Xmat = model.matrix(~Year, dat.northern)
coefs = data.frame(mod.northern) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
dat.northern = cbind(dat.northern,
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
)
save(dat.northern, file='data/modelled/dat.northern.RData')
save(dat.all.northern, file='data/modelled/dat.all.northern.RData')
save(mod.northern, file='data/modelled/mod.northern.RData')
load(file='data/modelled/dat.northern.RData')
load(file='data/modelled/dat.all.northern.RData')
load(file='data/modelled/mod.northern.RData')

### Calculate annual change
dat.northern = data.frame(Location='Northern',Year=sort(unique(dat.all.northern$Year)), N=length(unique(dat.all.northern$REEF_NAME)))
Xmat = model.matrix(~Year, dat.northern)
coefs = data.frame(mod.northern) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
contr.sequen=multcomp::contrMat(table(dat.northern$Year), type='Sequen')
Fit=Fit %*% t(contr.sequen)
dat.northern.diff = cbind(data.frame(Location='Northern'),
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
) %>% dplyr::rename(Years=X1) %>% mutate(Year=stringr::str_sub(Years, start=1, end=4))
save(dat.northern.diff, file='data/modelled/dat.northern.diff.RData')

## Central
dat.all.central = dat.all %>% filter(Location=='Central') %>% droplevels
#dat.all.central$Cvr1 = as.integer(as.vector(dat.all.central$Cover) * dat.all.central$Tows)
#dat.all.central$Cvr0 = dat.all.central$Tows - dat.all.central$Cvr1
mod.central <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|REEF_NAME), data=dat.all.central, family=binomial,iter=5000,warmup=2500,chains=3,cores=3)

dat.central = data.frame(Location='Central',Year=unique(dat.all.central$Year), N=length(unique(dat.all.central$REEF_NAME))) 
Xmat = model.matrix(~Year, dat.central)
coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
dat.central = cbind(dat.central,
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
)
save(dat.central, file='data/modelled/dat.central.RData')
save(dat.all.central, file='data/modelled/dat.all.central.RData')
save(mod.central, file='data/modelled/mod.central.RData')
load(file='data/modelled/dat.central.RData')
load(file='data/modelled/dat.all.central.RData')
load(file='data/modelled/mod.central.RData')

### Calculate annual change
dat.central = data.frame(Location='Central',Year=sort(unique(dat.all.central$Year)), N=length(unique(dat.all.central$REEF_NAME)))
Xmat = model.matrix(~Year, dat.central)
coefs = data.frame(mod.central) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
contr.sequen=multcomp::contrMat(table(dat.central$Year), type='Sequen')
Fit=Fit %*% t(contr.sequen)
dat.central.diff = cbind(data.frame(Location='Central'),
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
) %>% dplyr::rename(Years=X1) %>% mutate(Year=stringr::str_sub(Years, start=1, end=4))
save(dat.central.diff, file='data/modelled/dat.central.diff.RData')


## Southern
dat.all.southern = dat.all %>% filter(Location=='Southern') %>% droplevels
#dat.all.southern$Cvr1 = as.integer(as.vector(dat.all.southern$Cover) * dat.all.southern$Tows)
#dat.all.southern$Cvr0 = dat.all.southern$Tows - dat.all.southern$Cvr1
mod.southern <-  stan_glmer(cbind(Cvr1,Cvr0) ~ Year+(Year|REEF_NAME), data=dat.all.southern, family=binomial,iter=5000,warmup=2500,chains=3,cores=3)

dat.southern = data.frame(Location='Southern',Year=unique(dat.all.southern$Year), N=length(unique(dat.all.southern$REEF_NAME))) 
Xmat = model.matrix(~Year, dat.southern)
coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
dat.southern = cbind(dat.southern,
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
)
save(dat.southern, file='data/modelled/dat.southern.RData')
save(dat.all.southern, file='data/modelled/dat.all.southern.RData')
save(mod.southern, file='data/modelled/mod.southern.RData')
load(file='data/modelled/dat.southern.RData')
load(file='data/modelled/dat.all.southern.RData')
load(file='data/modelled/mod.southern.RData')

### Calculate annual change
dat.southern = data.frame(Location='Southern',Year=sort(unique(dat.all.southern$Year)), N=length(unique(dat.all.southern$REEF_NAME)))
Xmat = model.matrix(~Year, dat.southern)
coefs = data.frame(mod.southern) %>% dplyr:::select(matches('^X.*|^Year.*'))
Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
contr.sequen=multcomp::contrMat(table(dat.southern$Year), type='Sequen')
Fit=Fit %*% t(contr.sequen)
dat.southern.diff = cbind(data.frame(Location='Southern'),
            plyr:::adply(Fit,2,function(x) {
    data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
})
) %>% dplyr::rename(Years=X1) %>% mutate(Year=stringr::str_sub(Years, start=1, end=4))
save(dat.southern.diff, file='data/modelled/dat.southern.diff.RData')


## ggplot(dat.southern, aes(y=mean, x=Year)) +
##     geom_ribbon(aes(ymin=lower, ymax=upper, x=as.numeric(Year)), fill='blue', alpha=0.2)+
##     geom_line(aes(x=as.numeric(Year))) +
##     theme_bw()
