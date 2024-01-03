## --------------------------------------------------------------------------------
## Script to calculate the geostrophic wind values for London/Paris/De Bilt
## 
## Written by Richard Cornes (ricorne@noc.ac.uk) on 9th June 2023
## 
## --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

setwd("~/projects/geowind2")
source("scripts/geowind_function.R")

get_data <- function(data.dir, stn, fill_days){
    ## Read the data from data.dir, QC values and merge in fill_days
    
    require(dplyr)
    require(lubridate)
    require(dataresqc)
    
    ifiles <- list.files(data_dir,stn,full.names = TRUE)


    get_qc <- function(str,val) as.numeric(strsplit(grep(paste0(val,"="),strsplit(str,"\\|")[[1]],value=TRUE),"=")[[1]][2])
    
    data <- lapply(ifiles,function(F){
        meta <- read_meta(F,parameter = c("lon","lat"))
        print(F)

        dat <- read_sef(F,all=TRUE) %>%
            mutate(lon=meta["lon"],lat=meta["lat"],stn=stn,qc=get_qc(Meta,"qc"),Value=if_else(qc>0,NA_real_,Value))%>%
            select(Year:Minute,Value,lon,lat,stn)

        dat
    }) %>%
    bind_rows() %>%
    mutate(date=make_datetime(Year,Month,Day,Hour,Minute))

    if(stn=="Paris") data <- mutate(data,date=round(date,units="hour"))

    data <- data %>%
        mutate(Year=year(date),Month=month(date),Day=day(date),Hour=hour(date),Minute=minute(date)) %>%
        full_join(fill_days) %>%
        mutate(date=make_datetime(Year,Month,Day,Hour,Minute)) %>%
        select(date,Value,lon,lat,stn) %>%
        arrange(date)
    
    data    
}

fillgaps <- function(idata, station, maxgap=2){
    ## Inteprolate missing values using a linear interpolation
    
    require(zoo)
    require(dplyr)

    data.zoo <- zoo(idata$Value,idata$date)
    data.zoo <- na.approx(data.zoo, maxgap = maxgap)
    ##data.zoo <- na.spline(data.zoo, maxgap = maxgap,method="natural")

    as_tibble(data.zoo) %>%
        mutate(time=index(data.zoo),station=station)
}

## Setup a list of hourly values
tmin <- 1747
tmax <- 2023
fill.days <- data.frame(date=seq(make_datetime(tmin,1,1,12,0),
                             make_datetime(tmax,12,31,12,0), by="hour")) %>%
    mutate(Year=year(date),Month=month(date),Day=day(date),
           Hour=hour(date),Minute=minute(date)) %>%
    select(-date) %>%
    mutate_all(as.integer)

## Read the data
data_dir <- "~/projects/geowindWEB/data/SEF"
london <- get_data(data_dir,"London",fill.days)
paris <- get_data(data_dir,"Paris",fill.days)
debilt <- get_data(data_dir,"DeBilt",fill.days)


## Remove duplicates (two cases in the London series)
london <- london[-which(duplicated(london$date)),]
debilt <- debilt[-which(duplicated(debilt$date)),]

## Get value nearest 9am (+/- 2 hours)
london <- london %>%
    mutate(Year=year(date),Month=month(date),Day=day(date),Minute=minute(date),
           Hour=hour(date),ref=abs(Hour-9)) %>%
    filter(ref<=2&!is.na(Value)) %>%
    group_by(Year,Month,Day) %>%
    slice(which.min(ref)) %>%
    ungroup() %>%
    mutate_at(vars(Year,Month,Day,Hour,Minute),as.integer)%>%
    full_join(fill.days) %>%
    mutate(date=make_datetime(Year,Month,Day,Hour,Minute)) %>%
    select(date,Value,lon,lat,stn) %>%
    arrange(date)
    
## Interpolate subdaily values to the missing hours
idat <- list(london=london,paris=paris,debilt=debilt)
idat <- map2(idat, names(idat), function(x,y) fillgaps(x,y,maxgap=24)) %>%
    bind_rows() %>%
    spread(station,value)

## Calculate geostrophic wind
UV <- lapply(1:nrow(idat), function(i)
    geowind.UV(c(idat$paris[i],idat$london[i],idat$debilt[i]),
               c(48.8167,51.4787,52.11),
               c(2.333,-0.44904,5.183)))

idat$U <- unlist(lapply(UV,function(x)x[1]))
idat$V <- unlist(lapply(UV,function(x)x[2]))
idat$c1 <- unlist(lapply(UV,function(x)x[3]))
idat$c2 <- unlist(lapply(UV,function(x)x[4]))
idat$c3 <- unlist(lapply(UV,function(x)x[5]))

idat <- idat %>%
    mutate(geowind=geowind.speed(U,V),geowind.dir=geowind.dir(U,V),
           geowindQC=if_else(geowind>50,2,0))

saveRDS(idat,"data/geodat_QC_nonoon.RDS")

