#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Figs. S8 and S9
## Plots of seasonal and annual medians of the U and V wind components
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_geowind.R
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(latticeExtra)

setwd("~/projects/geowind2")

yrs <- c(1748,2023)

geodat <- readRDS("data/geodat_QC.RDS") %>%
    filter(minute(time)==0,hour(time)==12,between(year(time),yrs[1],yrs[2])) %>%
    mutate(U=if_else(!geowindQC==0,NA_real_,U),
           V=if_else(!geowindQC==0,NA_real_,V))

vector.seas <- geodat %>%
    mutate(month=month(time),year=year(time),
                      season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")), 
           year=ifelse(month%in%12,year+1,year)) %>%
    group_by(season,year) %>%
    summarise(U=median(U,na.rm = TRUE),V=median(V,na.rm = TRUE), n=length(which(!is.na(geowind)))/n()) %>%
    gather(component,data,-season,-year,-n) %>%
    mutate(data=if_else(n>0.8,data,NA_real_),component=factor(component),data.abs=abs(data)) %>%
    select(-n) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring"))) %>%
    na.omit()

plot.seas <- xyplot(data~year|season,data=vector.seas,type = "l",group=component,
                    auto.key = list(title="Component",columns=1,lines=TRUE,points=FALSE,cex=0.6,x=0.04,y=0.43),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    par.settings = list(superpose.line = list(pch = 19, cex = 1),
                                        par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm"))),
                    panel = function(x,y,...){
                        panel.refline(h=0, col="black",lwd=2)
                        panel.xyplot(x,y,...)}
                    )+
     glayer(panel.smoother(y ~ x+I(x^2), method = "lm",...))

pdf("plots/UV_season.pdf",width=7,height=6)
print(plot.seas)
dev.off()

vector.annual <- geodat %>%
    mutate(year=year(time)) %>%
    group_by(year) %>%
    summarise(U=median(U,na.rm = TRUE),V=median(V,na.rm = TRUE), n=length(which(!is.na(geowind)))/n()) %>%
    gather(component,data,-year,-n) %>%
    mutate(data=if_else(n>0.8,data,NA_real_),component=factor(component),data.abs=abs(data)) %>%
    select(-n) %>%
    na.omit()

plot.year <- xyplot(data~year,data=vector.annual,type = "l",group=component,
                    auto.key = list(title="Component",columns=1,lines=TRUE,points=FALSE,cex=0.6,x=0.8,y=0.1),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    par.settings = list(superpose.line = list(pch = 19, cex = 1),
                                        par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm"))),
                    panel = function(x,y,...){
                        panel.refline(h=0, col="black",lwd=2)
                        panel.xyplot(x,y,...)}
                    )+
     glayer(panel.smoother(y ~ x+I(x^2), method = "lm",...))

pdf("plots/UV_year.pdf",width=7,height=6)
print(plot.year)
dev.off()
