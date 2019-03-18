## ---- prelim
library(tidyverse)
library(lsmeans)
library(glmmTMB)
library(scales) # for brewer_pal
source('CoralTrends_functions.R')
## ----
pre.start=1985
pre.end = 2011

YearCuts = c(1985,2001,2011,2016)
CutLabels = c('Distant Past', 'Past', 'Recent')
CutLabels = c('1986-2001', '2002-2011', '2012-2016')
## Start with table of disturbances
## ---- LoadData
load(file='data/modelled/bleaching.full_3Zone.RData')
load(file='data/modelled/cots.full_3Zone.RData')
load(file='data/modelled/cyclones.full_3Zone.RData')
load(file='data/processed/all.reefs.cyclones.RData')
## ----
## New analysis

dist.table = list()

cyclones.full = all.reefs  %>% left_join(cyclones.full) %>% filter(!is.na(CYCLONEcat))
## spread the data so that for each reef/year there is a binary response for each category
cyclones.binary = cyclones.full %>%
    mutate(CYCLONEcat=as.factor(CYCLONEcat)) %>%
    bind_cols %>% data.frame(model.matrix(~-1+CYCLONEcat, data=.)) %>%
    mutate(CYCLONEany=ifelse(CYCLONEcat0==1, 0,1),
           CYCLONEsevere=ifelse((CYCLONEcat2+CYCLONEcat3)<1,0,1),
           time = REPORT_YEAR-min(REPORT_YEAR),
           Zone=factor(Zone, levels=c('Northern','Central','Southern')))

# ---- helperFunctions
calcFreqs = function(mod) {
    l1=lsmeans(mod, specs='time', by='Zone',type='response', at=list(time=0))
    list(Intercepts=l1)
    ##ideally, it would be good to be able to nominate the family to use here
    ##since the backtransform for slopes would just be exp rather than the
    ##inverse of a logit.  Add 1 to the slopes and intervals....
    ##SLOPES ARE ON A log odd-ratio scale.  WE SHOULD exp them so that they represent
    ## the factor change per year (or if them multiply by 100, the percent change per year)
    lt=lstrends(mod, specs='Zone', var='time')
    l2=test(lt)
    list(Intercepts=as.data.frame(summary(l1)) %>% full_join(test(l1)) %>% mutate(p.value=round(p.value,3)),
         slopes=as.data.frame(summary(lt)) %>% full_join(l2) %>% mutate(p.value=round(p.value,3)))
}

calcFreqs.matrix = function(mod) {
    dat=recover.data.glmmTMB(mod)#mod$frame
    form=formula(delete.response(terms(mod)))
    ## Slopes - actually rates (change in probability of being impacted per year)
    newdata = data.frame(time=1, Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
    Xmat = model.matrix(form, data=newdata)
    Xmat[,c(1,3,4)]=0
    coefs = fixef(mod)[[1]]
    (fit=as.vector(coefs %*% t(Xmat)))
    SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
    q=qnorm(0.975) #asymptotic (z test)
    l1=data.frame(fit=exp(fit), lower=exp(fit-q*SE), upper=exp(fit+q*SE))
    
    ## Intercepts - probabilty of being impacted at time 0 (1985)
    newdata = data.frame(time=0, Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
    Xmat = model.matrix(form, data=newdata)
    coefs = fixef(mod)[[1]]
    fit=as.vector(coefs %*% t(Xmat))
    SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
    q=qnorm(0.975) #asymptotic (z test)
    l2=data.frame(fit=binomial()$linkinv(fit), lower=binomial()$linkinv(fit-q*SE), upper=binomial()$linkinv(fit+q*SE))

    list(Intercept=l2, Slope=l1)
}

ACF.glmmTMB <- 
    function (object, maxLag, resType = c("pearson", "response", 
                                          "deviance","raw"), re=names(object$modelInfo$reTrms$cond$flist[1]),...) 
{
    resType <- match.arg(resType)
    res <- resid(object, type = resType)
    res = split(res,object$modelInfo$reTrms$cond$flist[[re]])
    if (missing(maxLag)) {
        maxLag <- min(c(maxL <- max(lengths(res)) - 1, as.integer(10 * 
                                                                  log10(maxL + 1))))
    }
    val <- lapply(res, function(el, maxLag) {
        N <- maxLag + 1L
        tt <- double(N)
        nn <- integer(N)
        N <- min(c(N, n <- length(el)))
        nn[1:N] <- n + 1L - 1:N
        for (i in 1:N) {
            tt[i] <- sum(el[1:(n - i + 1)] * el[i:n])
        }
        array(c(tt, nn), c(length(tt), 2))
    }, maxLag = maxLag)
    val0 <- rowSums(sapply(val, function(x) x[, 2]))
    val1 <- rowSums(sapply(val, function(x) x[, 1]))/val0
    val2 <- val1/val1[1L]
    z <- data.frame(lag = 0:maxLag, ACF = val2)
    attr(z, "n.used") <- val0
    class(z) <- c("ACF", "data.frame")
    z
}


recover.data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    recover.data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}

lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                               mode = "asymptotic", component="cond", ...) {
    if (mode != "asymptotic") stop("only asymptotic mode is available")
    if (component != "cond") stop("only tested for conditional component")
    if (missing(vcov.)) 
        V <- as.matrix(vcov(object)[[component]])
    else V <- as.matrix(.my.vcov(object, vcov.))
    dfargs = misc = list()
    if (!is.null(object$modelInfo$family)) {
        fam = object$modelInfo$family$family
        misc$tran = object$modelInfo$family$link
        misc$inv.lbl = "response"
        if (!is.na(pmatch(fam, "binomial"))) 
            misc$inv.lbl = "prob"
        else if (!is.na(pmatch(fam, "poisson"))) 
            misc$inv.lbl = "rate"
    }
    #misc = lsmeans:::.std.link.labels(object$modelInfo$family, misc)
    if (mode == "asymptotic") {
        dffun = function(k, dfargs) NA
    }
    ## use this? misc = .std.link.labels(family(object), misc)
    contrasts = attr(model.matrix(object), "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)[[component]]
    if (length(bhat) < ncol(X)) {
        kept = match(names(bhat), dimnames(X)[[2]])
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)[[component]]
        modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }
    else nbasis = estimability::all.estble
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc,...)
}


plotEffects = function(mod, firstYear,individualReefs=NULL, ribbonFillColor='blue') {
    dat=recover.data.glmmTMB(mod)#mod$frame
    tt <- terms(mod)
    Terms <- delete.response(tt)

    newdata = NULL
    for (z in unique(dat$Zone)) {
        dat1 = dat %>% filter(Zone==z)
        pts = with(dat1, expand.grid(Zone=z,
                                     time=seq(min(time), max(time), len=100)))
        newdata = rbind(newdata,pts)
    }
    newdata = newdata %>% mutate(Zone=factor(Zone, levels=c('Northern','Central','Southern')))
    #newdata = expand.grid(time=seq(min(dat$time), max(dat$time), len=100), Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
    m <- model.frame(Terms, newdata)
    Xmat <- model.matrix(Terms, m)
    #coefs = fixef(mod)[[1]]
    #fit=as.vector(coefs %*% t(X))
    
    #dat=recover.data.glmmTMB(mod)#mod$frame
    #form=formula(delete.response(terms(mod)))
    #newdata = expand.grid(time=seq(min(dat$time), max(dat$time), len=100), Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
    #Xmat = model.matrix(form, data=newdata)
    coefs = fixef(mod)[[1]]
    fit=as.vector(coefs %*% t(Xmat))
    SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
    q=qnorm(0.975) #asymptotic (z test)
    newdata = cbind(newdata, data.frame(fit=binomial()$linkinv(fit), lower=binomial()$linkinv(fit-q*SE), upper=binomial()$linkinv(fit+q*SE))) %>%
        mutate(Date=firstYear + time)

    if (!is.null(individualReefs)) {
        individualReefs = individualReefs %>% mutate(Date=firstYear + time, Time=factor(time))
        
        g1=ggplot(newdata, aes(y=fit, x=Date)) +
            geom_blank()+
            geom_point(data=individualReefs, aes(y=fit), color='grey') +
            geom_boxplot(data=individualReefs, aes(y=fit, group=Time), outlier.shape=NA) +
            facet_grid(~Zone) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill=ribbonFillColor, alpha=0.5, color=NA) +
            geom_line(color=ribbonFillColor) +
            scale_x_continuous('', expand=c(0,0), limits=c(1984,2020))+
            scale_y_continuous('Pr(impact)', limits=c(0,1.00))+
            theme_classic()+
            theme(strip.background = element_blank(), plot.margin=unit(c(0,2,0,1),'lines'),
                  panel.spacing=unit(1,'lines'))
    } else {
        g1=ggplot(newdata, aes(y=fit, x=Date)) +
            geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color=NA) +
            geom_line() +
            facet_grid(~Zone) +
            scale_y_continuous('Pr(impact)')+
            theme_classic()
    }
    g1
}
modelEachReef = function(mod,form=CYCLONEany~time, dat=mod$frame) {
    df = list()
    for (i in 1:length(unique(dat$REEF_NAME))) {
        d1=dat %>% filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
        if (nrow(d1)>3) {
            newmod = glm(form,data=d1, family=binomial)
            d1=data.frame(REEF_NAME=unique(d1$REEF_NAME), Zone=unique(d1$Zone),time=seq(min(d1$time),max(d1$time),by=1))
            df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
        }
    }
    df = do.call('rbind',df)
}
modelEachReefAR1 = function(mod,form=CYCLONEany~time) {
    dat=recover.data.glmmTMB(mod) %>% full_join(mod$frame)
    wch=grep('poly',colnames(dat))
    dat=dat[,-wch]
    df = list()
    for (i in 1:length(unique(dat$REEF_NAME))) {
        d1=dat %>% dplyr::filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
        newmod = glmmTMB(form,data=d1, family=binomial)
        d1=data.frame(REEF_NAME=unique(d1$REEF_NAME), Zone=unique(d1$Zone),time=seq(min(d1$time),max(d1$time),by=1))
        df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
    }
    df = do.call('rbind',df)
}


## ----


# Start off by exploring very short-term changes with splines.
library(mgcv)
cyclones.any.gamm = gamm(CYCLONEany ~ s(REPORT_YEAR, by=Zone), random=list(REEF_NAME=~1),
                           data=cyclones.binary, family=binomial())
summary(cyclones.any.gamm$gam)
plot(cyclones.any.gamm$gam, pages=1, shift=fixef(cyclones.any.gamm$lme)[[1]], trans=binomial()$linkinv, ylim=c(-5,5))

## GAM's over fit.  We really want to be able to say whether the long-term frequency of
## various disturbances has changed over time - GAMS are too short term.


library(glmmTMB)
cyclones.any.glmmTMB = glmmTMB(CYCLONEany ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
#cyclones.any.glmmTMB = glmmTMB(CYCLONEany ~ poly(time,3)*Zone + (1|REEF_NAME),
#                               data=cyclones.binary, family=binomial())
cyclones.any.glmmTMB1 = glmmTMB(CYCLONEany ~ time*Zone + (Zone|REEF_NAME),
                                data=cyclones.binary, family=binomial())
anova(cyclones.any.glmmTMB,cyclones.any.glmmTMB1)
## random intercept models fine..

plot(ACF(cyclones.any.glmmTMB, resType="pearson"), alpha=0.05)
## no evidence of temporal autocorrelation)

#cyclones.any.glmmTMB2 = glmmTMB(CYCLONEany ~ time*Zone + (1|REEF_NAME) + ar1(-1+time|Zone/REEF_NAME),
#                                data=cyclones.binary, family=binomial())

save(cyclones.any.glmmTMB, file='data/modelled/cyclones.any.glmmTMB.RData')

## ---- LoadCycloneModel
load(file='data/modelled/cyclones.any.glmmTMB.RData')
## ----
summary(cyclones.any.glmmTMB)
## Note, in the following, I would prefer the slopes where 1+slope estimate
## The underlying lstrends function backtransforms from logit (to produce odds) rather than
## log (to produce odds ratio).
## Ratio would be more intuitive on the natural scale as we can then say that
## for every one unit increase in time, the probability of experiencing an event
## increases by a factor of ...
## When backtransformed to odds, it does not have an interpretation.
## ---- CyclonesAny
calcFreqs(cyclones.any.glmmTMB)
calcFreqs.matrix(cyclones.any.glmmTMB)
df=modelEachReef(cyclones.any.glmmTMB,form=CYCLONEany~time)
d1=plotEffects(cyclones.any.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='output/figures/Disturbances_Cyclones.any.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.any.png',d1,width=7, height=3, dpi=300)
## ----


## ---- CyclonesSevere - This is the one used..
dist.table[['Cyclones']] = list()
cyclones.severe.glmmTMB = glmmTMB(CYCLONEsevere ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
plot(ACF(cyclones.severe.glmmTMB, resType="pearson"), alpha=0.05)
calcFreqs(cyclones.severe.glmmTMB)
calcFreqs.matrix(cyclones.severe.glmmTMB)
dist.table[['Cyclones']][['Intercept']] = calcFreqs(cyclones.severe.glmmTMB)[[1]] %>% dplyr::select(-time) %>% mutate(Disturbance='Cyclones', Stat='Intercept') %>% dplyr::select(Disturbance,Stat,everything())
dist.table[['Cyclones']][['Slope']] = calcFreqs(cyclones.severe.glmmTMB)[[2]] %>% mutate(Disturbance='Cyclones', Stat='Slope') %>% dplyr::select(Disturbance,Stat,everything())                                

df=modelEachReef(cyclones.severe.glmmTMB,form=CYCLONEsevere~time)
d1=plotEffects(cyclones.severe.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df, ribbonFillColor=brewer_pal(palette='Blues')(4)[4])
g.severe.cyclones=d1
ggsave(file='output/figures/Disturbances_Cyclones.severe.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.severe.png',d1,width=7, height=3, dpi=300)
## ----
save(cyclones.severe.glmmTMB, file='data/modelled/cyclones.severe.glmmTMB.RData')



ggplot(cyclones.binary, aes(y=CYCLONEcat3, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
cyclones.cat3.glmmTMB = glmmTMB(CYCLONEcat3 ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
plot(ACF(cyclones.cat3.glmmTMB, resType="pearson"), alpha=0.05)
## ---- CyclonesCat3
calcFreqs(cyclones.cat3.glmmTMB)
calcFreqs.matrix(cyclones.cat3.glmmTMB)
df=modelEachReef(cyclones.cat3.glmmTMB,form=CYCLONEcat3~time)
d1=plotEffects(cyclones.cat3.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='output/figures/Disturbances_Cyclones.cat3.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.cat3.png',d1,width=7, height=3, dpi=300)
## ----
save(cyclones.cat3.glmmTMB, file='data/modelled/cyclones.cat3.glmmTMB.RData')


ggplot(cyclones.binary, aes(y=CYCLONEcat2, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
cyclones.cat2.glmmTMB = glmmTMB(CYCLONEcat2 ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
plot(ACF(cyclones.cat2.glmmTMB, resType="pearson"), alpha=0.05)
## ---- CyclonesCat2
calcFreqs(cyclones.cat2.glmmTMB)
calcFreqs.matrix(cyclones.cat2.glmmTMB)
df=modelEachReef(cyclones.cat2.glmmTMB,form=CYCLONEcat2~time)
d1=plotEffects(cyclones.cat2.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='output/figures/Disturbances_Cyclones.cat2.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.cat2.png',d1,width=7, height=3, dpi=300)
## ----
save(cyclones.cat2.glmmTMB, file='data/modelled/cyclones.cat2.glmmTMB.RData')


##polynomials
cyclones.any.glmmTMB3 = glmmTMB(CYCLONEany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cyclones.binary, family=binomial())
save(cyclones.any.glmmTMB3, file='data/modelled/cyclones.any.glmmTMB3.RData')

plot(ACF(cyclones.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cyclones.any.glmmTMB3)
## ---- CyclonesAnyPoly
df=modelEachReef(cyclones.any.glmmTMB3,form=CYCLONEany~poly(time,3),dat=cyclones.binary)
d1=plotEffects(cyclones.any.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='output/figures/Disturbances_Cyclones.any_polygons.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.any_polygons.png',d1,width=7, height=3, dpi=300)
## ----

cyclones.severe.glmmTMB3 = glmmTMB(CYCLONEsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                   data=cyclones.binary, family=binomial())
save(cyclones.severe.glmmTMB3, file='data/modelled/cyclones.severe.glmmTMB3.RData')
plot(ACF(cyclones.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cyclones.severe.glmmTMB3)
## ---- CyclonesSeverePoly
df=modelEachReef(cyclones.severe.glmmTMB3,form=CYCLONEsevere~poly(time,3),dat=cyclones.binary)
d1=plotEffects(cyclones.severe.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='output/figures/Disturbances_Cyclones.severe_polygons.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Cyclones.severe_polygons.png',d1,width=7, height=3, dpi=300)
## ----




## COTS
cots.full = all.reefs  %>% left_join(cots.full %>% distinct) %>% filter(!is.na(COTScat))
## spread the data so that for each reef/year there is a binary response for each category
##For some reason, some reefs (e.g. '16017S') have two different COTScat in a particular year (1994)
## To correct for this, I will give them the max category.
cots.binary = cots.full %>%
    mutate(COTScat=factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    group_by(REEF_NAME,REEF_ID,Zone,Latitude,Longitude,REPORT_YEAR) %>%
    summarize(COTScat=levels(COTScat)[max(as.numeric(COTScat))]) %>%
    ungroup %>%
    mutate(COTScat=factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    bind_cols %>% data.frame(model.matrix(~-1+COTScat, data=.)) %>%
    mutate(COTSany=ifelse(COTScatZero==1, 0,1),
           COTSsevere=ifelse((COTScatIO+COTScatAO)<1,0,1),
           time = REPORT_YEAR-min(REPORT_YEAR),
           Time=as.factor(time),
           Zone=factor(Zone, levels=c('Northern','Central','Southern')))
save(cots.binary, file='data/modelled/cots.binary.RData')

g1 = ggplot(cots.binary, aes(y=COTSany, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/cots_any.pdf', g1,width=10, height=3)

#cots.any.glmmTMB = glmmTMB(COTSany ~ time*Zone + (1|REEF_NAME) + ar1(as.factor(time)-1|REEF_NAME),
#                           data=cots.binary, family=binomial())
cots.any.glmmTMB = glmmTMB(COTSany ~ time*Zone + (1|REEF_NAME)+ ar1(-1+factor(time)|REEF_NAME),
                           data=cots.binary, family=binomial())
summary(cots.any.glmmTMB)
plot(ACF(cots.any.glmmTMB, resType="pearson"), alpha=0.05)
## ---- CotsAny
calcFreqs(cots.any.glmmTMB)
calcFreqs.matrix(cots.any.glmmTMB)
## ----
#df=modelEachReefAR1(cots.any.glmmTMB,form=COTSany~time+ar1(-1+factor(time)|REEF_NAME))
df=modelEachReef(cots.any.glmmTMB,form=COTSany~time,dat=cots.binary)
d1=plotEffects(cots.any.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='output/figures/Disturbances_COTS_any.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_any.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.any.glmmTMB.pdf', d1,width=7, height=3)
## ----
save(cots.any.glmmTMB, file='data/modelled/cots.any.glmmTMB.RData')

##polynomials
cots.any.glmmTMB3 = glmmTMB(COTSany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.any.glmmTMB3, file='data/modelled/cots.any.glmmTMB3.RData')

plot(ACF(cots.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.any.glmmTMB3)
## ---- CotsAnyPoly
#df=modelEachReefAR1(cots.any.glmmTMB3,form=COTSany~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.any.glmmTMB3.df.RData')
df=modelEachReef(cots.any.glmmTMB3,form=COTSany~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.any.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
                                        #ggsave(file='figures/cots.any.glmmTMB3.pdf', g1,width=10, height=3)
ggsave(file='output/figures/Disturbances_COTS_any_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_any_poly.png',d1,width=7, height=3, dpi=300)
## ----



g1=ggplot(cots.binary, aes(y=COTSsevere, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/cots_severe.pdf', g1,width=10, height=3)

#cots.severe.glmmTMB = glmmTMB(COTSsevere ~ time*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
#                              data=cots.binary, family=binomial())
cots.severe.glmmTMB = glmmTMB(COTSsevere ~ time*Zone + (1|REEF_NAME),
                               data=cots.binary, family=binomial())
plot(ACF(cots.severe.glmmTMB, resType="pearson"), alpha=0.05)
## ---- CotsSevere
calcFreqs(cots.severe.glmmTMB)
calcFreqs.matrix(cots.severe.glmmTMB)
## ----
#df=modelEachReefAR1(cots.severe.glmmTMB,form=COTSsevere~time+ar1(-1+factor(time)|REEF_NAME))
df=modelEachReef(cots.severe.glmmTMB,form=COTSsevere~time, dat=cots.binary)
d1=plotEffects(cots.severe.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
g.severe.cots = d1
ggsave(file='output/figures/Disturbances_COTS_severe.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_severe.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.severe.glmmTMB.pdf', g1,width=10, height=3)
## ----
save(cots.severe.glmmTMB, file='data/modelled/cots.severe.glmmTMB.RData')

##polynomials
#cots.severe.glmmTMB3 = glmmTMB(COTSsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
#                               data=cots.binary, family=binomial())
cots.severe.glmmTMB3 = glmmTMB(COTSsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.severe.glmmTMB3, file='data/modelled/cots.severe.glmmTMB3.RData')

plot(ACF(cots.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.severe.glmmTMB3)
## ---- CotsSeverePoly
#df=modelEachReefAR1(cots.severe.glmmTMB3,form=COTSsevere~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.severe.glmmTMB3.df.RData')
df=modelEachReef(cots.severe.glmmTMB3,form=COTSsevere~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.severe.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='output/figures/Disturbances_COTS_severe_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_severe_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.severe.glmmTMB3.pdf', g1,width=10, height=3)
## ----



## CotsAO - This is the one used..
dist.table[['COTS']] = list()
g1=ggplot(cots.binary, aes(y=COTScatAO, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/cots_catAO.pdf', g1,width=10, height=3)

cots.catAO.glmmTMB = glmmTMB(COTScatAO ~ time*Zone + (1|REEF_NAME),# + ar1(-1+time|REEF_NAME),
                             data=cots.binary, family=binomial())
save(cots.catAO.glmmTMB, file='data/modelled/cots.catAO.glmmTMB.RData')
plot(ACF(cots.catAO.glmmTMB, resType="pearson"), alpha=0.05)
## ---- CotsCatAO
calcFreqs(cots.catAO.glmmTMB)
calcFreqs.matrix(cots.catAO.glmmTMB)
dist.table[['COTS']][['Intercept']] = calcFreqs(cots.catAO.glmmTMB)[[1]] %>% dplyr::select(-time) %>% mutate(Disturbance='COTS', Stat='Intercept') %>% dplyr::select(Disturbance,Stat,everything())
dist.table[['COTS']][['Slope']] = calcFreqs(cots.catAO.glmmTMB)[[2]] %>% mutate(Disturbance='COTS', Stat='Slope') %>% dplyr::select(Disturbance,Stat,everything())

## ----
#df=modelEachReefAR1(cots.catAO.glmmTMB,form=COTScatAO~time+ar1(-1+time|REEF_NAME))
df=modelEachReef(cots.catAO.glmmTMB,form=COTScatAO~time)
d1=plotEffects(cots.catAO.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
g.AO.cots = d1
ggsave(file='output/figures/Disturbances_COTS_AO.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_AO.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.catAO.glmmTMB.pdf', g1,width=10, height=3)
## ----

##polynomials
cots.catAO.glmmTMB3 = glmmTMB(COTScatAO ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.catAO.glmmTMB3, file='data/modelled/cots.catAO.glmmTMB3.RData')

plot(ACF(cots.catAO.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.catAO.glmmTMB3)
## ---- CotsCatAOPoly
#df=modelEachReefAR1(cots.catAO.glmmTMB3,form=COTScatAO~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.catAO.glmmTMB3.df.RData')
df=modelEachReef(cots.catAO.glmmTMB3,form=COTScatAO~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.catAO.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='output/figures/Disturbances_COTS_AO_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_COTS_AO_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.catAO.glmmTMB3.pdf', g1,width=10, height=3)
## ----


## ggplot(cots.binary, aes(y=COTScatIO, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
## cots.catIO.glmmTMB = glmmTMB(COTScatIO ~ time*Zone + (1|REEF_NAME) + ar1(-1+time|REEF_NAME),
##                              data=cots.binary, family=binomial())
## cots.catIO.glmmPQL = glmmPQL(COTScatIO ~ time*Zone, random=~1|REEF_NAME,
##                              data=cots.binary, family=binomial(), correlation=corAR1(form=~time|REEF_NAME))
## cots.catIO.glmmPQL = glmmPQL(COTScatIO ~ poly(time,3)*Zone, random=~1|REEF_NAME,
##                              data=cots.binary, family=binomial(), correlation=corAR1(form=~time|REEF_NAME))
## plot(allEffects(cots.catIO.glmmPQL))
## plot(ACF(cots.catIO.glmmTMB, resType="pearson"), alpha=0.05)
## ## ---- CotsIO
## calcFreqs(cots.catIO.glmmTMB)
## calcFreqs.matrix(cots.catIO.glmmTMB)
## df=modelEachReefAR1(cots.catIO.glmmTMB,form=COTScatIO~poly(time,3)+ar1(-1+time|REEF_NAME))
## plotEffects(cots.catIO.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df)
## ## ----
## save(cots.catIO.glmmTMB, file='data/modelled/cots.catIO.glmmTMB.RData')



## ##polynomials
## cyclones.any.glmmTMB3 = glmmTMB(CYCLONEany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
##                                 data=cyclones.binary, family=binomial())
## save(cyclones.any.glmmTMB3, file='data/modelled/cyclones.any.glmmTMB3.RData')

## plot(ACF(cyclones.any.glmmTMB3, resType="pearson"), alpha=0.05)
## summary(cyclones.any.glmmTMB3)
## ## ---- CyclonesAnyPoly
## df=modelEachReef(cyclones.any.glmmTMB3,form=CYCLONEany~poly(time,3),dat=cyclones.binary)
## plotEffects(cyclones.any.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
## ## ----

## cyclones.severe.glmmTMB3 = glmmTMB(CYCLONEsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
##                                    data=cyclones.binary, family=binomial())
## save(cyclones.severe.glmmTMB3, file='data/modelled/cyclones.severe.glmmTMB3.RData')
## plot(ACF(cyclones.severe.glmmTMB3, resType="pearson"), alpha=0.05)
## summary(cyclones.severe.glmmTMB3)
## ## ---- CyclonesSeverePoly
## df=modelEachReef(cyclones.severe.glmmTMB3,form=CYCLONEsevere~poly(time,3),dat=cyclones.binary)
## plotEffects(cyclones.severe.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
## ## ----




## bleaching
## mike supplied some bleaching data to fill in the gaps between aerial surveys
load(file='data/modelled/bleaching.merge_3Zone.RData')
#bleaching.full = all.reefs  %>% left_join(bleaching.merge %>% distinct) %>% filter(!is.na(BLEACHINGcat))
bleaching.full = bleaching.merge %>% filter(!is.na(Zone)) %>% droplevels %>%
    group_by(REEF_ID) %>% arrange(REPORT_YEAR) %>% ungroup %>%
    mutate(Zone=factor(Zone, levels=c('Northern','Central','Southern')))

bleaching.binary = bleaching.full %>% filter(!is.na(BLEACHINGcat)) %>%
    mutate(BLEACHINGcat=as.factor(BLEACHINGcat)) %>%
    bind_cols %>% data.frame(model.matrix(~-1+BLEACHINGcat, data=.)) %>%
    mutate(BLEACHINGany=ifelse(BLEACHINGcat0==1, 0,1),
           BLEACHINGsevere=ifelse((BLEACHINGcat3+BLEACHINGcat4)<1,0,1),
           time = REPORT_YEAR-min(REPORT_YEAR),
           Zone=factor(Zone, levels=c('Northern','Central','Southern')))
save(bleaching.binary, file='data/modelled/bleaching.binary.RData')

g1=ggplot(bleaching.binary, aes(y=BLEACHINGany, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/bleaching_any.pdf', g1,width=10, height=3)

bleaching.any.glmmTMB = glmmTMB(BLEACHINGany ~ time*Zone + (1|REEF_NAME),
                                data=bleaching.binary, family=binomial())
save(bleaching.any.glmmTMB, file='data/modelled/bleaching.any.glmmTMB.RData')

summary(bleaching.any.glmmTMB)
plot(ACF(bleaching.any.glmmTMB, resType="pearson"), alpha=0.05)
## ---- BleachingAny
calcFreqs(bleaching.any.glmmTMB)
calcFreqs.matrix(bleaching.any.glmmTMB)
## ----
df=modelEachReef(bleaching.any.glmmTMB,form=BLEACHINGany~time) %>% mutate(Zone=factor(Zone, levels=c('Northern','Central','Southern')))
d1=plotEffects(bleaching.any.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
ggsave(file='output/figures/Disturbances_Bleaching_any.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Bleaching_any.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/bleaching.any.glmmTMB.pdf', g1,width=10, height=3)
## ----


##polynomials
bleaching.any.glmmTMB3 = glmmTMB(BLEACHINGany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=bleaching.binary, family=binomial())
save(bleaching.any.glmmTMB3, file='data/modelled/bleaching.any.glmmTMB3.RData')

plot(ACF(bleaching.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(bleaching.any.glmmTMB3)
## ---- BleachingAnyPoly
#df=modelEachReefAR1(bleaching.any.glmmTMB3,form=BLEACHINGany~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/bleaching.any.glmmTMB3.df.RData')
df=modelEachReef(bleaching.any.glmmTMB3,form=BLEACHINGany~poly(time,3),dat=bleaching.binary)
d1=plotEffects(bleaching.any.glmmTMB3, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
ggsave(file='output/figures/Disturbances_Bleaching_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Bleaching_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/bleaching.any.glmmTMB3.pdf', g1,width=10, height=3)
## ----




g1=ggplot(bleaching.binary, aes(y=BLEACHINGsevere, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/bleaching_severe.pdf', g1,width=10, height=3)

## ---- Bleaching severe - This is the one used..
dist.table[['Bleaching']] = list()
bleaching.severe.glmmTMB = glmmTMB(BLEACHINGsevere ~ time*Zone + (1|REEF_NAME),
                                   data=bleaching.binary, family=binomial())
save(bleaching.severe.glmmTMB, file='data/modelled/bleaching.severe.glmmTMB.RData')

summary(bleaching.severe.glmmTMB)
plot(ACF(bleaching.severe.glmmTMB, resType="pearson"), alpha=0.05)
## ---- BleachingSevere
calcFreqs(bleaching.severe.glmmTMB)
calcFreqs.matrix(bleaching.severe.glmmTMB)
dist.table[['Bleaching']][['Intercept']] = calcFreqs(bleaching.severe.glmmTMB)[[1]] %>% dplyr::select(-time) %>% mutate(Disturbance='Bleaching', Stat='Intercept') %>% dplyr::select(Disturbance,Stat,everything())
dist.table[['Bleaching']][['Slope']] = calcFreqs(bleaching.severe.glmmTMB)[[2]] %>% mutate(Disturbance='Bleaching', Stat='Slope') %>% dplyr::select(Disturbance,Stat,everything())  
## ----
#df=modelEachReefAR1(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time+ar1(-1+time|REEF_NAME))
df=modelEachReef(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time)
#df=modelEachReef(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time)
d1=plotEffects(bleaching.severe.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
g.severe.bleaching = d1
ggsave(file='output/figures/Disturbances_Bleaching_severe.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Bleaching_severe.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/bleaching.severe.glmmTMB.pdf', g1,width=10, height=3)
## ----

##polynomials - does not converge..
#bleaching.severe.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
#                               data=bleaching.binary, family=binomial())
bleaching.severe.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,2,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=bleaching.binary, family=binomial())
save(bleaching.severe.glmmTMB3, file='data/modelled/bleaching.severe.glmmTMB3.RData')

plot(ACF(bleaching.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(bleaching.severe.glmmTMB3)
## ---- BleachingSeverePoly
#df=modelEachReefAR1(bleaching.severe.glmmTMB3,form=BLEACHINGsevere~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/bleaching.severe.glmmTMB3.df.RData')
df=modelEachReef(bleaching.severe.glmmTMB3,form=BLEACHINGsevere~poly(time,3),dat=bleaching.binary)
d1=plotEffects(bleaching.severe.glmmTMB3, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
ggsave(file='output/figures/Disturbances_Bleaching_severe_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_Bleaching_severe_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/bleaching.severe.glmmTMB3.pdf', g1,width=10, height=3)
## ----


## Save the disturbance tables
save(dist.table, file='data/modelled/dist.table.RData')
write.csv(do.call('rbind',lapply(dist.table, `[[`, 'Intercept')), file='data/modelled/dist.table.intercepts.csv', quote=FALSE, row.names=FALSE)
write.csv(do.call('rbind',lapply(dist.table, `[[`, 'Slope')), file='data/modelled/dist.table.slopes.csv', quote=FALSE, row.names=FALSE)

## Combine all disturbances
disturb.binary = cyclones.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,CYCLONEany,CYCLONEsevere) %>%
    full_join(cots.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,COTSany, COTSsevere)) %>%
    full_join(bleaching.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,BLEACHINGany, BLEACHINGsevere) %>% mutate(REEF_ID=as.numeric(as.character(REEF_ID)))) %>%
    distinct

disturb.binary = disturb.binary %>% filter(!is.na(CYCLONEany), !is.na(BLEACHINGany), !is.na(COTSany)) %>%
    mutate(DISTURBany=ifelse((CYCLONEany+COTSany+BLEACHINGany)>0,1,0),
           DISTURBsevere = ifelse((CYCLONEsevere+COTSsevere+BLEACHINGsevere)>0,1,0),
           time = REPORT_YEAR-min(REPORT_YEAR),
           Zone=factor(Zone, levels=c('Northern','Central','Southern')),
           REEF_NAME=REEF_ID
           )
save(disturb.binary, file='data/modelled/disturb.binary.RData')
g1=ggplot(disturb.binary, aes(y=DISTURBany, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/disturb_any.pdf', g1,width=10, height=3)

disturb.any.glmmTMB = glmmTMB(DISTURBany ~ time*Zone + (1|REEF_NAME),
                              data=disturb.binary, family=binomial())
save(disturb.any.glmmTMB, file='data/modelled/disturb.any.glmmTMB.RData')

summary(disturb.any.glmmTMB)
plot(ACF(disturb.any.glmmTMB, resType="pearson"), alpha=0.05)
## ---- DisturbAny
calcFreqs(disturb.any.glmmTMB)
calcFreqs.matrix(disturb.any.glmmTMB)
## ----
df=modelEachReef(disturb.any.glmmTMB,form=DISTURBany~time)
d1=plotEffects(disturb.any.glmmTMB, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='output/figures/Disturbances_any.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_any.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.any.glmmTMB.pdf', g1,width=10, height=3)
## ----

##polynomials
disturb.any.glmmTMB3 = glmmTMB(DISTURBany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=disturb.binary, family=binomial())
save(disturb.any.glmmTMB3, file='data/modelled/disturb.any.glmmTMB3.RData')
plot(ACF(disturb.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(disturb.any.glmmTMB3)
## ---- DisturbAnyPoly
df=modelEachReef(disturb.any.glmmTMB3,form=DISTURBany~poly(time,3),dat=disturb.binary)
d1=plotEffects(disturb.any.glmmTMB3, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='output/figures/Disturbances_any_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_any_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.any.glmmTMB3.pdf', g1,width=10, height=3)
## ----


## Severe disturbances
#disturb.binary = cyclones.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,CYCLONEsevere) %>%
#    full_join(cots.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,COTSsevere)) %>%
#    full_join(bleaching.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,BLEACHINGsevere) %>% mutate(REEF_ID=as.numeric(as.character(REEF_ID)))) %>%
#    distinct

#disturb.binary = disturb.binary %>% filter(!is.na(CYCLONEsevere), !is.na(BLEACHINGsevere), !is.na(COTSsevere)) %>%
#    mutate(DISTURBsevere=ifelse((CYCLONEsevere+COTSsevere+BLEACHINGsevere)>0,1,0),
#           time = REPORT_YEAR-min(REPORT_YEAR),
#           Zone=factor(Zone, levels=c('Northern','Central','Southern')),
#           REEF_NAME=REEF_ID
#           )

g1=ggplot(disturb.binary, aes(y=DISTURBsevere, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='figures/disturb_severe.pdf', g1,width=10, height=3)

disturb.severe.glmmTMB = glmmTMB(DISTURBsevere ~ time*Zone + (1|REEF_NAME),
                                 data=disturb.binary, family=binomial())
save(disturb.severe.glmmTMB, file='data/modelled/disturb.severe.glmmTMB.RData')
summary(disturb.severe.glmmTMB)
plot(ACF(disturb.severe.glmmTMB, resType="pearson"), alpha=0.05)
## ---- DisturbSevere
calcFreqs(disturb.severe.glmmTMB)
calcFreqs.matrix(disturb.severe.glmmTMB)
## ----
df=modelEachReef(disturb.severe.glmmTMB,form=DISTURBsevere~time)
d1=plotEffects(disturb.severe.glmmTMB, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
g.severe.any = d1
ggsave(file='output/figures/Disturbances_severe.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_severe.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.severe.glmmTMB.pdf', g1,width=10, height=3)
## ----

##polynomials
disturb.severe.glmmTMB3 = glmmTMB(DISTURBsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=disturb.binary, family=binomial())
save(disturb.severe.glmmTMB3, file='data/modelled/disturb.severe.glmmTMB3.RData')
plot(ACF(disturb.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(disturb.severe.glmmTMB3)
## ---- DisturbSeverePoly
df=modelEachReef(disturb.severe.glmmTMB3,form=DISTURBsevere~poly(time,3),dat=disturb.binary)
d1=plotEffects(disturb.severe.glmmTMB3, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='output/figures/Disturbances_severe_poly.pdf',d1,width=7, height=3)
ggsave(file='output/figures/Disturbances_severe_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.severe.glmmTMB3.pdf', g1,width=10, height=3)
## ----


## Now form a single multipanel figure with the following
## 1. All severe disturbances
## 2. COTS severe disturbances (IO and AO)
## 3. Cyclone severe disturbances (Cat 2,3)
## 4. Bleaching severe disturbances (Cat 3,4)
gg=grid.arrange(g.severe.any + ggtitle('All severe disturbances'),
             g.severe.cots + ggtitle('A. cf. solaris IO and AO'),
             g.severe.cyclones + ggtitle('Cyclones Hs cat 2 and 3'),
             g.severe.bleaching + ggtitle('Bleaching cat 3 and 4'), ncol=1)

ggsave(file='output/figures/Disturbances_severe_compilation.pdf',grid.draw(gg),width=9, height=12)
ggsave(file='output/figures/Disturbances_severe_compilation.png',grid.draw(gg),width=9, height=12, dpi=300)

## Now form a single multipanel figure with the following
## 1. COTS severe disturbances (IO and AO)
## 2. Cyclone severe disturbances (Cat 2,3)
## 3. Bleaching severe disturbances (Cat 3,4)
gg=grid.arrange(
             g.severe.cots + ggtitle('A. cf. solaris IO and AO'),
             g.severe.cyclones + ggtitle('Cyclones Hs cat 2 and 3'),
             g.severe.bleaching + ggtitle('Bleaching cat 3 and 4'), ncol=1)

ggsave(file='output/figures/Disturbances_severe_compilation_noAll.pdf',grid.draw(gg),width=9, height=9)
ggsave(file='output/figures/Disturbances_severe_compilation_noAll.png',grid.draw(gg),width=9, height=9, dpi=300)

## ----


## Now combine with disturbances from CoralTrend_trend_manta_3Zone.R
load(file='data/processed/disturbanceBars-gg.RData')
g.severe.any1=g.severe.any +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    #scale_y_continuous(expression(Pr(impact)),lim=c(0,1.00)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,0),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
gg2 = ggplot_gtable(ggplot_build(g.severe.any1))
#gg3=gg2

#gg2=gg3
gg$widths[4] =  gg2$widths[4]
#gg2$widths=gg$widths
gg2$widths[6] =  gg$widths[7]
gg2$widths[8] =  gg$widths[7] 
grid.arrange(gg,gg2)
ggsave(file='output/figures/Disturbances_severe_full_compilation.pdf',grid.arrange(gg,gg2),width=9, height=6)
ggsave(file='output/figures/Disturbances_severe_full_compilation.png',grid.arrange(gg,gg2),width=9, height=6, dpi=300)


##Finally, a figure that combines the disturbance bars from CoralTrends_trend_manta_3Zone.R with the severe version of each major disturbance.
load(file='data/processed/disturbanceBars-gg.RData')
g.severe.cots1= g.severe.cots + ggtitle('c) A. cf. solaris IO and AO') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))

g.severe.cots1= g.AO.cots + ggtitle('c) A. cf. solaris AO') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))


g.severe.cyclones1= g.severe.cyclones + ggtitle('d) Cyclones Hs cat 2 and 3') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cyclones2 = ggplot_gtable(ggplot_build(g.severe.cyclones1))

g.severe.bleaching1= g.severe.bleaching + ggtitle('b) Bleaching cat 3 and 4') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1985,2020))+
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
                                        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.bleaching2 = ggplot_gtable(ggplot_build(g.severe.bleaching1))
gg$widths[4] =  g.severe.bleaching2$widths[4]
g.severe.bleaching2$widths[6] =  gg$widths[7]
g.severe.bleaching2$widths[8] =  gg$widths[7]
g.severe.cots2$widths[6] =  gg$widths[7]
g.severe.cots2$widths[8] =  gg$widths[7]
g.severe.cyclones2$widths[6] =  gg$widths[7]
g.severe.cyclones2$widths[8] =  gg$widths[7]

gg1=grid.arrange(
    g.severe.bleaching2,
    g.severe.cots2,
    g.severe.cyclones2,
    ncol=1,
    left=textGrob(expression(Probability~of~being~impacted~phantom("(")), vjust=.375,rot=90, gp=gpar(fontsize=16)))
    #left=textGrob(expression(Probability~of~being~impacted), rot=90, gp=gpar(fontsize=16)))
    #left=expression(Probability~of~being~impacted~phantom("(")))
    #left='Probability of being impacted')
grid.arrange(gg,gg1, heights=c(1,3))
ggsave(file='output/figures/Disturbances_severe_compilation.pdf',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8)
ggsave(file='output/figures/Disturbances_severe_compilation.png',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8, dpi=300)


## Zip up some stuff for Mike
zip('output/DisturbanceFrequencyFigures4Mike.zip', list.files('output/figures', pattern='Disturbances.*.pdf', full.names=TRUE))

load(file='data/processed/disturbanceBars-gg.RData')
g.severe.any1=g.severe.any +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    #scale_y_continuous(expression(Pr(impact)),lim=c(0,1.00)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,0),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
gg2 = ggplot_gtable(ggplot_build(g.severe.any1))
    






