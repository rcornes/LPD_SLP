#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce the Geowind trend tables
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
## --------------------------------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(zyp)
library(broom)
library(xtable)

setwd("~/projects/geowind2")

pc_file <- "data/geowind_seasonal_percentiles.csv"
DF <- read.csv(pc_file) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")),quantile=factor(pc)) %>%
    na.omit()


idates <- list(c(1871,2023),c(1748,2023),c(1950,2023))
cnt <- 1

for (rng in idates){

    results <- DF %>%
    filter(between(year,rng[1],rng[2])) %>%
    select(year,season,quantile,Q) %>%
    group_by(season,quantile) %>%
    nest() %>%
    mutate(model=map(data,~zyp.zhang(.x$Q,.x$year)),
           tidied=map(model,tidy)) %>%
    unnest(tidied) %>%
    select(-data,-model) %>%
        filter(names%in%c("lbound","trend","ubound","sig")) %>%
        spread(names,x) %>%
        mutate_at(.vars=vars(lbound,trend,ubound),.funs=funs(.*100)) %>%
        mutate_at(.vars=vars(lbound,trend,ubound),.funs=funs(round(.,2))) %>%
        mutate(sig=round(sig,4)) %>%
        select(season,quantile,trend,lbound,ubound)

    if(cnt==1){
        caption <- sprintf("Trends (per century) and 95 percent confidence level (indicated by the lbound and ubound values) in seasonal geostrophic wind percentiles over the period %s--%s. ",rng[1],rng[2])
        label <- "tab:trend1"
    } else {
        caption <- sprintf("As Table \\ref{tab:trend1} except for the period %s-%s.",rng[1],rng[2])
        label <- paste0("tab:trend",cnt)
    }

    print(xtable(results,caption=caption,label=label),file=sprintf("tables/trends_table_%s--%s.tex",rng[1],rng[2]),include.rownames=FALSE)

    cnt <- cnt+1
}
