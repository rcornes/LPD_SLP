#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. 1
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles.R, make_boot.R
## --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(latticeExtra)
library(splines)
library(tactile)

setwd("~/projects/geowind2")

yrs <- c(1747,2023)

## Read percentiles
pc_file <- "data/geowind_seasonal_percentiles.csv"
boot_file <- "data/geowind_seas_boot_noon.csv"

DF <- read.csv(pc_file) %>%
    filter(pc%in%c("5th","50th","95th"),between(year,yrs[1],yrs[2])) %>%
    mutate(pc=factor(pc,levels=c("5th","50th","95th"))) %>%
    na.omit()

DF.boot <- read.csv(boot_file) %>%
    filter(ci%in%c("0.05","0.95")) %>%
    gather(pc,boot,-year,-season,-ci) %>%
    mutate(pc=factor(pc,levels=paste0("G",c(0.05,0.5,0.95)),labels=c("5th","50th","95th")),
           ci=factor(ci,levels=c("0.05","0.95"),labels=c("lower.ci","upper.ci")),
           boot=round(boot,3)) %>%
    spread(ci,boot)

DF <- left_join(DF,DF.boot) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")))

## Plot the percentiles
my3cols <- c("#FC4E07", "#2E9FDF","#E7B800")

DF1 <- filter(DF,between(year,1747,1774))
DF2 <- filter(DF,between(year,1774,2023))
xlim <- c(1740,2030)

plot.seas <- xyplot(Q~year|season,data=DF1,type = "l",group=pc, prepanel = prepanel.ci,
                    auto.key = list(title="Percentile",columns=1,lines=TRUE,points=FALSE,
                                    cex=0.6,x=0.25,y=0.45),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    lower=DF1$lower.ci,upper=DF1$upper.ci,
                    yscale.components = yscale.components.subticks,
                    lwd=.5,xlim=xlim,ylim=c(0,38),
                    panel=function(...,col=my3cols){
                        panel.ci(...,alpha=0.25,grid=TRUE)
                        panel.xyplot(...,alpha=0.5)
                    })

plot.seas2 <- xyplot(Q~year|season,data=DF2,type = "l",group=pc, prepanel = prepanel.ci,
                    lower=DF2$lower.ci, upper=DF2$upper.ci, xlim=xlim,ylim=c(0,38),
                    lwd=.5,yscale.components = yscale.components.subticks,
                    panel=function(...){
                        panel.ci(...,alpha=0.3)
                        panel.xyplot(...)
                    })+
    glayer(panel.smoother(y ~  ns(x,5), method = "lm",se=FALSE,lower=NULL,upper=NULL,...))

## Save to file
pdf("plots/fig_2.pdf",width=6.5,height=6)
plot(plot.seas+plot.seas2)
dev.off()
