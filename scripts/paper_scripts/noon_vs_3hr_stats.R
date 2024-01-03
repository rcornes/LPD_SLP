#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce 3-hourly vs. noon geowind stats
## Tables of RMSE, trend and linear regression
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-06
##
## Prerequisites scripts: make_percentiles.R, make_percentiles_3hr.R
## --------------------------------------------------------------------------------


library(tidyverse)
library(broom)
library(lubridate)
library(zyp)
library(xtable)

setwd("~/projects/geowind2")

## Read percentiles
DF.noon <- read.csv("data/geowind_seasonal_percentiles.csv") %>%
    filter(pc%in%c("5th","50th","95th"),year>=1950) %>%
    mutate(series="Noon") %>%
    na.omit()

DF.3hr <- read.csv("data/geowind_seasonal_percentiles_3hr.csv") %>%
    filter(pc%in%c("5th","50th","95th")) %>%
    mutate(series="3 hourly") %>%
    na.omit()

DF <- rbind(DF.noon,DF.3hr) %>%
    mutate(season=factor(season,levels=c("Autumn","Summer","Spring","Winter")),
           pc=factor(pc,levels=c("5th","50th","95th")))

## Trends
DF.trends <- DF %>%
    group_by(season,pc,series) %>%
    nest() %>%
    mutate(model=map(data,~zyp.zhang(.x$Q,.x$year)),
           tidied=map(model,tidy)) %>%
    unnest(tidied) %>%
    select(-data,-model) %>%
        filter(names%in%c("lbound","trend","ubound")) %>%
        spread(names,x) %>%
        mutate_at(.vars=vars(lbound,trend,ubound),.funs=funs(.*100)) %>%
    mutate_at(.vars=vars(lbound,trend,ubound),.funs=funs(round(.,2))) %>%
    mutate(trend_rng=sprintf("%s (%s - %s)",trend,lbound,ubound)) %>%
    select(trend_rng) %>%
    spread(series,trend_rng) %>%
    as.data.frame()

caption <- "Linear trends and 95pctile uncertainty range of the seasonal L-P-Dtri time series over the period 1950-2021 using noon-only and three--hourly sampling."



print(xtable(DF.trends,caption=caption), file = "tables/trends_3hr_noon_table.tex",include.rownames = FALSE)

DF.regress <- DF %>%
    select(season,year,pc,Q,series) %>%
    spread(series,Q) %>%
    group_by(season,pc) %>%
    nest() %>%
    mutate(model=map(data,~lm(.x$Noon~.x$`3 hourly`)),
           tidied=map(model,glance)) %>%
    unnest(tidied) %>%
    select(season,pc,adj.r.squared) %>%
    mutate(adj.r.squared=round(adj.r.squared,2)) %>%
    spread(season,adj.r.squared)

caption <- "Linear regression of 3-hourly and noon-only percentiles"

print(xtable(DF.regress,caption=caption), file = "tables/regress_3hr_noon_table.tex",include.rownames = FALSE)

rmse <- function(x,y) sqrt(mean((x-y)^2,na.rm=TRUE))

DF.rmse <- DF %>%
    select(season,year,pc,Q,series) %>%
    spread(series,Q) %>%
    group_by(season,pc) %>%
    nest() %>%
    mutate(rmse=map(data,~rmse(.x$Noon, .x$`3 hourly`))) %>%
    unnest(rmse) %>%
    select(season,pc,rmse) %>%
    mutate(rmse=round(rmse,2)) %>%
    spread(season,rmse)

caption <- "Root, mean-squared error of the noon-only L-P-Dtri season quantiles against the three-hourly sampling series."
print(xtable(DF.rmse,caption=caption), file = "tables/rmse_3hr_noon_table.tex",include.rownames = FALSE)
