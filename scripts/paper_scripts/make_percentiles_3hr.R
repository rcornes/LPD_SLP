#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce the 3-hourly geowind seasonal percentiles
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
## --------------------------------------------------------------------------------


library(dplyr)
library(tidyr)
library(lubridate)

setwd("~/projects/geowind2")

yrs <- c(1950,2023) #Year range

## Select only noon values and QC values
geodat <- readRDS("data/geodat_QC.RDS") %>%
    filter(between(year(time),yrs[1],yrs[2]),minute(time)==0,hour(time)%in%seq(0,23,3)) %>%
    mutate(geowind=if_else(!geowindQC==0,NA_real_,geowind))

## Calculate Seasonal Percentiles
geowind.seas <- geodat %>%
    mutate(month=month(time),year=year(time),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")), 
           year=ifelse(month%in%12,year+1,year))%>%
    group_by(season,year)%>%
    reframe(Q=quantile(geowind,c(0.05,0.5,0.90,0.95,0.99),na.rm = TRUE),pc=c("5th","50th","90th","95th","99th"),n=length(which(!is.na(geowind)))/n())%>%
    mutate(Q=if_else(n<0.8,NA_real_,Q),n=round(n,2),Q=round(Q,3)) %>%
    arrange(year,season)


## Write to file
write.csv(geowind.seas,"data/geowind_seasonal_percentiles_3hr.csv",row.names=FALSE,quote=FALSE)
