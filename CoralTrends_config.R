########################################################################
## This module loads a cofiguration file that contains paramaters for ##
## settings to be applied across all aspects of the analyses          ##
## This file is a text file with key=value pairs                      ##
##                                                                    ##
## Specifically, the pairs are:                                       ##
## Size:                    the number of bootstrapp samples          ##
## Timeseries.plot.width    the width of a timeseries plot            ##
## EndDate                  the last day (date) of the current focal  ##
##                            year (must be in YYYY-MM-DD format)     ##
## Timeseries.plot.height   the height of a timeseries plot           ##
## FocalYear                the current report card year              ##
## StartDate                the lower date range cutoff (must be in   ##
##                            YYY-MM-DD format)                       ##
########################################################################


CoralTrends_tryCatch(
    {
        files <- list.files('data')
        if(!any(grepl('^primary$',files))) system('mkdir data/primary')
        if(!any(grepl('^processed$',files))) system('mkdir data/processed')

        files <- list.files()
        if(!any(grepl('^figures$',files))) system('mkdir figures')

        files <- list.files()
        if(!any(grepl('^output$',files))) system('mkdir output')
        if(!any(grepl('^output$',files))) system('cd output; mkdir figures;')

        files <- list.files()
        if(!any(grepl('^logs$',files))) system('mkdir logs')
        return=NULL
    }, 'logs/all.log','--Config--',msg='configure necessary folders', return=NULL)

CoralTrends_tryCatch(
    {
        config = readLines('parameters/CoralTrends.conf')
        config = gsub('(.*Date)=(.*)','\\1=as.Date(\'\\2\')',config)
        eval(parse(text=config))
        return=NULL
    }, 'logs/all.log','--Config--',msg='load general configurations', return=NULL)


