#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Supplementary Figure 6
## Plot of Annual Percentiles of Geostrophic Wind
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-11
##
## Prerequisites scripts: make_percentiles.R, make_boot.R
## --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(latticeExtra)
library(tactile)
library(splines)

setwd("~/projects/geowind2")

## Read percentiles
pc_file <- "data/geowind_annual_percentiles.csv"
boot_file <- "data/geowind_annual_boot_noon.csv"

DF <- read.csv(pc_file) %>%
    filter(pc%in%c("5th","50th","95th")) %>%
    mutate(pc=factor(pc,levels=c("5th","50th","95th"))) %>%
    na.omit()

DF.boot <- read.csv(boot_file) %>%
    filter(ci%in%c("0.05","0.95")) %>%
    gather(pc,boot,-year,-ci) %>%
    mutate(pc=factor(pc,levels=paste0("G",c(0.05,0.5,0.95)),labels=c("5th","50th","95th")),
           ci=factor(ci,levels=c("0.05","0.95"),labels=c("lower.ci","upper.ci"))) %>%
    spread(ci,boot)

DF <- left_join(DF,DF.boot)

## Plot the percentiles
my3cols <- c("#FC4E07", "#2E9FDF","#E7B800")

DF1 <- filter(DF,between(year,1747,1774))
DF2 <- filter(DF,between(year,1774,2023))
xlim <- c(1740,2030)

plot.ann <- xyplot(Q~year,data=DF1,type = "l",group=pc, prepanel = prepanel.ci,
                   auto.key = list(title="Percentile",columns=3,lines=TRUE,points=FALSE,cex=0.6,x=0.2,y=1),
                    xlab="Year",ylab=expression(paste("Geostrophic Wind Speed (",m, " ", s^-1,")", sep="")),
                    lower=DF1$lower.ci,upper=DF1$upper.ci,
                    lwd=.5,xlim=xlim,
                    panel=function(...,col=my3cols){
                        panel.ci(...,alpha=0.1,grid=TRUE)
                        panel.xyplot(...,alpha=0.5)
                    })

plot.ann2 <- xyplot(Q~year,data=DF2,type = "l",group=pc, prepanel = prepanel.ci,
                    lower=DF2$lower.ci,upper=DF2$upper.ci,xlim=xlim,
                    lwd=.5,
                    panel=function(...){
                        panel.ci(...,alpha=0.25)
                        panel.xyplot(...)
                    })+
     glayer(panel.smoother(y ~  ns(x,5), method = "lm",lower=NULL,upper=NULL,...))


TS.plot <- plot.ann+plot.ann2

## Save to file
pdf("plots/annual_percentiles.pdf")
print(TS.plot)
dev.off()
