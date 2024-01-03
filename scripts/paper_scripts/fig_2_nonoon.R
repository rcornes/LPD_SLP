#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. 1
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles_nonoon.R
## --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(latticeExtra)
library(gridExtra)
library(splines)

setwd("~/projects/geowind2")

yrs <- c(1747,2023)

## Read percentiles
pc_file <- "data/geowind_seasonal_percentiles_nonoon.csv"

DF <- read.csv(pc_file) %>%
    filter(pc%in%c("5th","50th","95th"),between(year,yrs[1],yrs[2])) %>%
    mutate(pc=factor(pc,levels=c("5th","50th","95th")),
           season=factor(season,levels=c("Summer","Autumn","Winter","Spring"))) %>%
    na.omit()

DF1 <- filter(DF,between(year,1747,1774))
DF2 <- filter(DF,between(year,1774,2023))
xlim <- c(1740,2030)

plot.seas <- xyplot(Q~year|season,data=DF1,type = "l",group=pc, #prepanel = prepanel.ci,
                    auto.key = list(title="Percentile",columns=1,lines=TRUE,points=FALSE,
                                    cex=0.6,x=0.25,y=0.45),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    lwd=.5,xlim=xlim,ylim=c(0,38),
                    yscale.components = yscale.components.subticks,
                    panel=function(...,col=my3cols){
                        ##panel.ci(...,alpha=0.1,grid=TRUE)
                        panel.xyplot(...,alpha=0.5)
                    })

plot.seas2 <- xyplot(Q~year|season,data=DF2,type = "l",group=pc, #prepanel = prepanel.ci,
                    lower=DF2$lower.ci, ylim=c(0,38),
                    lwd=.5,,yscale.components = yscale.components.subticks,
                    panel=function(...){
                        ##panel.ci(...,alpha=0.25)
                        panel.xyplot(...)
                    })+
     glayer(panel.smoother(y ~  ns(x,5), method = "lm",se=FALSE,lower=NULL,upper=NULL,...))


TS.plot <- plot.seas+plot.seas2

## Save to file
pdf("plots/fig_2_nonoon.pdf",width=6,height=6)
print(TS.plot)
dev.off()
