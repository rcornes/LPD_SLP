#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. S3
## Plot of percentiles 1950-2023 in percentiles from noon and 3-hourly data
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-06
##
## Prerequisites scripts: make_percentiles.R, make_percentiles_3hr.R
## --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(tactile)
library(latticeExtra)
library(splines)

setwd("~/projects/geowind2")

## Read percentiles
DF.noon <- read.csv("data/geowind_seasonal_percentiles.csv") %>%
    filter(pc%in%c("5th","50th","95th"),year>=1950) %>%
    mutate(series="Noon") %>%
    na.omit()

DF.noon.boot <- read.csv("data/geowind_seas_boot_noon.csv") %>%
    filter(ci%in%c("0.05","0.95")) %>%
    gather(pc,boot,-year,-season,-ci) %>%
    mutate(pc=factor(pc,levels=paste0("G",c(0.05,0.5,0.95)),labels=c("5th","50th","95th")),
           ci=factor(ci,levels=c("0.05","0.95"),labels=c("lower.ci","upper.ci"))) %>%
    spread(ci,boot)

DF.noon <- left_join(DF.noon,DF.noon.boot) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")))

DF.3hr <- read.csv("data/geowind_seasonal_percentiles_3hr.csv") %>%
    filter(pc%in%c("5th","50th","95th")) %>%
    mutate(series="3 hourly") %>%
    na.omit()

DF.3hr.boot <- read.csv("data/geowind_seas_boot_3hr.csv") %>%
    filter(ci%in%c("0.05","0.95")) %>%
    gather(pc,boot,-year,-season,-ci) %>%
    mutate(pc=factor(pc,levels=paste0("G",c(0.05,0.5,0.95)),labels=c("5th","50th","95th")),
           ci=factor(ci,levels=c("0.05","0.95"),labels=c("lower.ci","upper.ci"))) %>%
    spread(ci,boot)

DF.3hr <- left_join(DF.3hr,DF.3hr.boot) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")))

DF <- rbind(DF.noon,DF.3hr) %>%
    mutate(season=factor(season,levels=c("Autumn","Summer","Spring","Winter")),
           pc=factor(pc,levels=c("5th","50th","95th")))

## Construct plot
plot.seas <- xyplot(Q~year|series+season,data=DF,type = "l",group=pc, prepanel = prepanel.ci,
                    auto.key = TRUE,
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    lower=DF$lower.ci,upper=DF$upper.ci,
                    lwd=.5,
                    par.settings = list(superpose.line = list(pch = 19, cex = 1),
                                        par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm"))),
                    panel=function(...,col=my3cols){
                        panel.ci(...,alpha=0.25,grid=TRUE)
                        panel.xyplot(...)
                    })+
     glayer(panel.smoother(y ~  ns(x,5), method = "lm",lower=NULL,upper=NULL,...))

## Save to file
pdf("plots/geowind_percentiles_noon_vs_3hr.pdf",height=10, width=7)
print(plot.seas)
dev.off()
