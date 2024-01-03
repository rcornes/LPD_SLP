#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. S5
## Plot of annual gale frequencies in geowind, UK Gales and Armagh series.
##
## This script also calculates the difference in freqs. in the geowind
## and Armagh documentary series over the periods 1800-20 and 1821-50
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles.R
##
## Other data sources:
## 1. The gale index is downloaded from the CRU website within the script
## 2. The Armargh series should be downloaded from:
##    http://climate.armagh.ac.uk/calibrated/storm/
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(latticeExtra)
library(splines)

setwd("~/projects/geowind2")

ann.freqs <- read_csv("data/geowind_annual_frequency.csv",col_select = c(year,N)) %>%
    mutate(series="English Channel Storms",source="English Channel Storms") %>%
    rename(data=N)

gales.in <- read_csv(url("https://crudata.uea.ac.uk/cru/data/lwt/webdocs/gale_index_1200hrs_UK.csv"),col_names=TRUE)
#gales.in.extended <- read_csv("~/projects/LWTs/20CRv3_LWTs_output/20CRv3_1806_2015_12hrs_outfile_UK.csv")

gales.ann <- gales.in %>%
    filter(!`Gale Index exceedances`%in%c("gales")) %>%
    rename(year=Year) %>%
    select(-`Gale Index exceedances`) %>%
    group_by(year) %>%
    summarise(data=n()) %>%
    mutate(series="UK Gales",source="UK Gales") %>%
    arrange(year)

arm.ann <- read.fwf("data/armargh/annual-diary",c(5,2,99), skip=12, na="nd", col.names = c("year","data","dummy")) %>%
    select(-dummy) %>%
    mutate(series="Armagh",source="Armagh diary")

arm.ins.ann <- read.fwf("data/armargh/annual-instrument",c(5,2,99), skip=11, na="nd", col.names = c("year","data","dummy"))%>%
    select(-dummy) %>%
    mutate(series="Armagh",source="Armagh Instrumental",data=as.numeric(data))

gg.ann.data <- bind_rows(arm.ann,ann.freqs,gales.ann,arm.ins.ann) %>%
    mutate(data=ifelse(is.na(data),0,data),source=if_else(is.na(source),"Instrumental",source))

ann.plot <- xyplot(data~year|source,data=gg.ann.data,type = "s",group=source,
                      xlab="Year",ylab="Frequency",lwd=0.25,
                      par.settings = list(par.main.text = list(font = 1, just = "left", x = grid::unit(18, "mm")),
                                          superpose.line=list(lwd=0.5)),
                    panel=function(...){
                        panel.abline(h=seq(0,40,5),v=seq(1800,2000,50),col="lightgrey",lwd=0.5)
                        panel.xyplot(...)
                    })+
     glayer(panel.smoother(y ~ ns(x,5), method = "lm",...))
 
pdf("plots/gales_annual.pdf",height=4,width=6)
print(ann.plot)
dev.off()

## Compare periods
gg.ann.data %>%
    filter(between(year,1800,1850)) %>%
    mutate(period=case_when(between(year,1800,1820)~"1800-1820",
                            between(year,1821,1850)~"1821-1850")) %>%
    group_by(source,period) %>%
    summarise(av=mean(data,na.rm=TRUE)) %>%
    spread(period,av) %>%
    mutate(diff=(`1800-1820`-`1821-1850`)/`1800-1820`*100)%>%
    filter(source%in%c("English Channel Storms","Armagh diary"))
