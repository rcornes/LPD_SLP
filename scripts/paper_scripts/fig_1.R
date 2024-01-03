#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce Fig. 1
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-12-30
##
## --------------------------------------------------------------------------------

library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(lubridate)
library(cowplot)
library(ggplot2)
library(metR)
library(tactile)

setwd("~/projects/geowind2")

## Plot ERA5 on stormiest two days
world <- ne_countries(scale = "medium", returnclass = "sf")

d1 <- raster("data/ERA5/era5_SLP_2501199012Z.nc")
d2 <- raster("data/ERA5/era5_SLP_1411199312Z.nc")
d1.date <- getZ(d1)
d2.date <- getZ(d2)
lon.lim <- c(-15,15.5)
lat.lim <- c(46,68)

d1 <- as.data.frame(d1,xy=TRUE) %>%
    as_tibble() %>%
    rename(d1=Mean.sea.level.pressure)

d2 <- as.data.frame(d2,xy=TRUE) %>%
    as_tibble() %>%
    rename(d2=Mean.sea.level.pressure)

dd <- inner_join(d1,d2) %>%
    gather(date,mslp,-x,-y) %>%
    mutate(date=factor(date,levels=c("d1","d2"),labels=c("a)","b)")),mslp=mslp/100)%>%
    filter(between(x,lon.lim[1],lon.lim[2]),between(y,lat.lim[1],lat.lim[2]))

map.storms <- ggplot()+
    geom_sf(data=world,linewidth=0.25)+
    geom_contour(aes(x,y,z=mslp),data=dd,breaks=seq(940,1030,4),linewidth=0.25)+
    theme_map(12)+
    coord_sf(xlim=lon.lim,ylim=lat.lim,expand=FALSE)+
    geom_segment(aes(x=2.333,xend=-0.44904,y=48.8167,yend=51.4787),data=world,colour="red",linewidth=0.1)+
    geom_segment(aes(x=-0.44904,xend=5.183,y=51.4787,yend=52.11),data=world,colour="red",linewidth=0.1)+
    geom_segment(aes(x=5.183,xend=2.333,y=52.109,yend=48.8167),data=world,colour="red",linewidth=0.1)+
    annotate("point", x = c(2.333,-0.44904,5.183), y = c(48.8167,51.4787,52.11), colour = "red", size = 1)+
    geom_text_contour(aes(x,y,z=mslp), data=dd,stroke = 0.2,size=2,breaks=seq(940,1030,8))+
    theme(strip.text.x=element_text(hjust = 0))+
    facet_wrap(~date)

pdf("plots/fig_1.pdf",width=5.5,height=3.5)
plot(map.storms)
dev.off()
