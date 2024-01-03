#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce bootstrap confidence estimates
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-10
##
## Prerequisites scripts: make_geowind.R
## --------------------------------------------------------------------------------

library(dplyr)
library(lubridate)
library(doParallel)

setwd("~/projects/geowind2/")

args <- commandArgs(trailingOnly=TRUE)

if (!length(args)==1) {
  stop("Provide time resolution (noon or 3hr)", call.=FALSE)
}

stopifnot(args[1]%in%c("noon","3hr"))

if(args[1]=="noon"){
    res <- 12
    year1 <- 1747
    year2 <- 2023
} else {
    res <- seq(0,23,3)
    year1 <- 1948
    year2 <- 2023
}

message(sprintf("Processing %s - %s", year1, year2))

## --------------------------------------------------------------------------------
## Hard-coded Parameters
pctiles <- c(0.05,0.5,0.95) #Quantiles to calculate
nsim <- 10000               #Number of bootstrap samples
nclus <- 4                  #Number of parallel cores to use
## --------------------------------------------------------------------------------

geowind.boot <- function(x, quant, B=10000, nfac=0.8){

    x <- x[!is.na(x)]
    n_obs <- length(x)
    nsamp <- n_obs*nfac
    ci <- c(0.05,0.5,0.95)

    if(all(is.na(x))){
        Mout <- matrix(NA_real_,nrow=length(ci),ncol=length(quant),dimnames = list(1:length(ci),paste0("G",quant)))
    } else {
        result_vec <- matrix(nrow=B,ncol=length(quant),dimnames = list(1:B,paste0("G",quant)))

        for(b in 1:B) {
            this_sample <- sample(x, size=nsamp, replace=TRUE)
            m <- quantile(this_sample,quant,na.rm=TRUE)
            result_vec[b,] <- m
        }

        Mout <- apply(result_vec,2,quantile,probs=ci)
    }
    
    Mout <- data.frame(Mout)
    Mout$ci <- ci

    Mout
}


geodat <- readRDS("data/geodat_QC.RDS")

## Read data and subset by years/times
geodat <- geodat %>%
    filter(minute(time)==0,hour(time)%in%res,between(year(time),year1,year2)) %>%
    mutate(geowind=if_else(!geowindQC==0,NA_real_,geowind))

geowind.seas <- geodat %>%
    mutate(month=month(time),year=year(time),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")), 
           year=ifelse(month%in%12,year+1,year))

## Run seasonal boostrap
cl <- makeCluster(nclus)
registerDoParallel(cl)

boot.seas <- foreach(Y=unique(geowind.seas$year),.combine="rbind")%:%
    foreach(seas=unique(geowind.seas$season),.combine="rbind",.packages="dplyr")%dopar%{
        dat <- filter(geowind.seas,year==Y,season==seas)
        odat <- geowind.boot(dat$geowind,pctiles,B=nsim)
        mutate(odat,year=Y,season=seas)
    }

write.csv(boot.seas, sprintf("data/geowind_seas_boot_%s.csv",args[1]),row.names=FALSE,quote=FALSE)

## Run annual boostrap
boot.year <- foreach(Y=unique(geowind.seas$year),.combine="rbind",.packages="dplyr")%dopar%{
        dat <- filter(geowind.seas,year==Y)
        odat <- geowind.boot(dat$geowind,pctiles,B=nsim)
        mutate(odat,year=Y)
    }

write.csv(boot.year, sprintf("data/geowind_annual_boot_%s.csv",args[1]),row.names=FALSE,quote=FALSE)
