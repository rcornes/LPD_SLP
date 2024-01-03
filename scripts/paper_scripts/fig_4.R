#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. 3
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_geowind.R
##
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(openair)

setwd("~/projects/geowind2")

geodat <- readRDS("data/geodat_QC.RDS") %>%
    rename(date=time) %>%
    filter(minute(date)==0,hour(date)==12) %>%
    mutate(geowind=if_else(!geowindQC==0,NA_real_,geowind),
           geowind.dir=if_else(!geowindQC==0,NA_real_,geowind.dir),
           month=month(date),year=year(date),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%c(3:5)~"Spring",
                            month%in%c(6:8)~"Summer",
                            month%in%c(9:11)~"Autumn"),
           year=if_else(month==12,year+1,year)) %>%
    mutate(period=case_when(between(year,1791,1820)~"1791-1820",
                            between(year,1990,2023)~"1990-2023"))

seas <- list("Winter","Spring","Summer","Autumn")
names(seas) <- seas #Inherit names in later lapplys
lims <- list("Winter"=seq(0,30,6),"Spring"=seq(0,20,4),"Summer"=seq(0,20,4),"Autumn"=seq(0,30,6))
labs <- list("Winter"="a)",
             "Spring"="b)",
             "Summer"="a) Summer (June - August) ",
             "Autumn"="b) Autumn (September - November)")

my.settings <- list(par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm")))

## Seasonal Plots
wr.plots <- lapply(seas, function(S){
    windRose(filter(geodat,season==S), ws = 'geowind', wd = 'geowind.dir',
             paddle = FALSE, border = FALSE, breaks = lims[[S]],
             key=list(height=1,space="bottom",
                      header=expression(paste("Geostrophic  Wind Speed (",m, " ", s^-1,")", sep="")),
                      footer=NULL,plot.style = c("ticks", "border")),
             type="period",
             par.settings=my.settings,
             cex=0.1,
             lwd=0.1,
             main=labs[[S]],
             angle=45,
             layout=c(2,1),
             annotate=FALSE,
             dig.lab=1,
             grid.line = 5)
})

pdf("plots/fig_4.pdf",width=6.5, height=9)
plot(wr.plots[["Winter"]],split=c(1,1,1,2))
plot(wr.plots[["Spring"]],split=c(1,2,1,2), newpage=FALSE)
dev.off()

## Annual plot
wr.plots.year <- windRose(geodat, ws = 'geowind', wd = 'geowind.dir',
                        paddle = FALSE, border = FALSE,
                        breaks = seq(0,30,6),
                        key=list(height=1,space="bottom",
                                 header=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                                 footer=NULL,plot.style = c("ticks", "border")),
                        type="period",
                        par.settings=my.settings,
                        cex=0.1,
                        lwd=0.1,
                        main="c) Annual",
                        angle=45,
                        layout=c(2,1),
                        annotate=FALSE,
                        dig.lab=1,
                        grid.line = 5)


pdf("plots/windrose_annual_summer_autumn.pdf",width=7, height=12)
plot(wr.plots[["Summer"]],split=c(1,1,1,3))
plot(wr.plots[["Autumn"]],split=c(1,2,1,3), newpage=FALSE)
plot(wr.plots.year,split=c(1,3,1,3), newpage=FALSE)
dev.off()

