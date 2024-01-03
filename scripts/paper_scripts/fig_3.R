#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. 2
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles.R
##
## Other data sources:
## 1. The gale index is downloaded from the CRU website within the script
## 2. The Armagh series should be downloaded from:
##    http://climate.armagh.ac.uk/calibrated/storm/
## --------------------------------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(broom)
library(gtools)
library(gridExtra)
library(latticeExtra)
library(splines)

setwd("~/projects/geowind2")

yrs <- c(1747,2023)

geodat.seas <- read_csv("data/geowind_seasonal_frequency.csv",col_select = c(year,season,N)) %>%
    filter(between(year,yrs[1],yrs[2])) %>%
    mutate(series="English Channel Storms",source="English Channel Storms") %>%
    rename(data=N)

gales.seas <- read_csv(url("https://crudata.uea.ac.uk/cru/data/lwt/webdocs/gale_index_1200hrs_UK.csv"),col_names=TRUE) %>%
    filter(`Gale Index exceedances`==c("severe gales")) %>%
    rename(year=Year) %>%
    mutate(season=case_when(Month%in%c(12,1,2)~"Winter",
                            Month%in%3:5~"Spring",
                            Month%in%6:8~"Summer",
                            Month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")),
           year=ifelse(Month%in%12,year+1,year)) %>%
    group_by(year,season) %>%
    summarise(data=n()) %>%
    mutate(series="UK Gales", source="UK Gales")

arm.seas <- read.fwf("data/armargh/seasonal-diary",c(4,rep(2,4)), skip=5,col.names=c("year","Winter","Spring","Summer","Autumn")) %>%
    gather(season,data,-year) %>%
    mutate(series="Armagh",source="Armagh Diary") %>%
    select(year,season,series,data,source)

arm.ins.seas <- read.fwf("data/armargh/seasonal-instrument",c(4,rep(2,4),100), skip=5,
                         col.names=c("year","Winter","Spring","Summer","Autumn","components")) %>%
    select(-components) %>%
    gather(season,data,-year) %>%
    mutate(series="Armagh",source="Armagh Instrumental") %>%
    select(year,season,series,data,source) %>%
    mutate(data=if_else(is.na(data),0,data))

gg.winter.data <- bind_rows(geodat.seas,gales.seas,arm.seas,arm.ins.seas) %>%
    mutate(data=ifelse(is.na(data),0,data),source=if_else(is.na(source),"Instrumental",source)) %>%
    filter(season=="Winter") %>%
    mutate(data=ifelse(is.na(data),0,data))

dat1 <- gg.winter.data %>%
    mutate(data=if_else(year>1774,NA_real_,data))
dat2 <- gg.winter.data %>%
    mutate(data=if_else(year<1774,NA_real_,data))

xlim <- c(1740,2030)
ylim <- c(0,20)

plot1 <- xyplot(data~year|source,data=dat1,type = "s",group=source,
                      xlab="Year",ylab="Frequency",lwd=0.25,xlim=xlim,ylim=ylim,
                      par.settings = list(par.main.text = list(font = 1, just = "left", x = grid::unit(18, "mm")),
                                          superpose.line=list(lwd=0.5)),
                    panel=function(...){
                        panel.abline(h=seq(0,20,5),v=seq(1800,2000,50),col="lightgrey",lwd=0.5)
                        panel.xyplot(...,alpha=0.5)
                    })

plot2 <- xyplot(data~year|source,data=dat2,type = "s",group=source,
                      xlab="Year",ylab="Frequency",lwd=0.25,xlim=xlim,ylim=ylim,
                      par.settings = list(par.main.text = list(font = 1, just = "left", x = grid::unit(18, "mm")),
                                          superpose.line=list(lwd=0.5)),
                    panel=function(...){
                        panel.abline(h=seq(0,20,5),v=seq(1800,2000,50),col="lightgrey",lwd=0.5)
                        panel.xyplot(...)
                    })+
    glayer(panel.smoother(y ~ ns(x,5), method = "lm",...))

TS.plot <- plot1+plot2
 
pdf("plots/fig_3.pdf",height=4,width=6)
print(TS.plot)
dev.off()

## Compare periods
gg.winter.data %>%
    filter(between(year,1800,1850)) %>%
    mutate(period=case_when(between(year,1800,1820)~"1800-1820",
                            between(year,1821,1850)~"1821-1850")) %>%
    group_by(source,period,season) %>%
    summarise(av=mean(data,na.rm=TRUE)) %>%
    spread(period,av) %>%
    mutate(diff=(`1800-1820`-`1821-1850`)/`1800-1820`*100)
