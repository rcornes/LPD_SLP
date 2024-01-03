#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Supplementary Figure 4
## Plot of 20CR-derived geostrophic wind values per season
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles.R
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(latticeExtra)

setwd("~/projects/geowind2")
source("scripts/geowind_function.R")

dat <- map(c("london","debilt","paris"), function(x){
    read.table(file.path("data",sprintf("%s_20CRv3.txt",x)),head=FALSE,col.names=c("year","month","day","time","value")) %>%
        as_tibble() %>%
        mutate(value=value/100) %>%
        rename(!!x:=value)
})

dat <- reduce(dat,inner_join)

UV <- lapply(1:nrow(dat), function(i)
    geowind.UV(c(dat$paris[i],dat$london[i],dat$debilt[i]),
                   c(48.857,51.507,52.109),
                   c(2.351,0.128,5.181)))

dat$U <- unlist(lapply(UV,function(x)x[1]))
dat$V <- unlist(lapply(UV,function(x)x[2]))

geodat <- dat %>%
    mutate(geowind=geowind.speed(U,V),geowind.dir=geowind.dir(U,V)) %>%
    mutate(season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Summer","Autumn","Winter","Spring")),
           year=ifelse(month%in%12,year+1,year))

## PLOTS
xmin <- 1806
xmax <- 2015

## Seasonal plot
geowind.seas <- geodat %>%
    mutate(time=as.character(time))%>%
    filter(time=="12:00:00") %>%
    group_by(year,season) %>%
        summarise(geowind.05=quantile(geowind,0.05,na.rm=TRUE),
                  geowind.5=quantile(geowind,0.5,na.rm=TRUE),
                  geowind.95=quantile(geowind,0.95,na.rm=TRUE),
                  N=length(which(is.na(geowind)))/length(geowind),
                  geowind.05=ifelse(N<0.2,geowind.05,NA),
                  geowind.5=ifelse(N<0.2,geowind.5,NA),
                  geowind.95=ifelse(N<0.2,geowind.95,NA))%>%
    dplyr::select(year,season,geowind.05,geowind.5,geowind.95)%>%
    gather(Percentile,data,-year,-season)%>%
    mutate(Percentile=factor(Percentile,levels=c("geowind.95","geowind.5","geowind.05"),
                             labels = c("95th","50th","5th"))) %>%
    dplyr::filter(between(year,xmin,xmax))

my3cols <- c("#FC4E07", "#2E9FDF","#E7B800")
plot.seas <- xyplot(data~year|season,data=geowind.seas,type = "l",group=Percentile,
                    auto.key = list(title="Percentile",columns=1,lines=TRUE,points=FALSE,cex=0.6,x=0.04,y=0.43),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    par.settings = list(superpose.line = list(pch = 19, cex = 1,col = my3cols),
                                        par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm"))))+
     glayer(panel.smoother(y ~ x+I(x^2), method = "lm",...))

pdf("plots/geowind_season_20CR_1806_2015.pdf",width=6,height=6)
print(plot.seas)
dev.off()



