#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce NAO/geowind correlation table
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_percentiles.R
## The NAO indices (apart from the EA index) are downloaded within the script.
## The EA index has been downloaded from: https://www.cpc.ncep.noaa.gov/data/teledoc/ea.shtml
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(broom)
library(xtable)

setwd("~/projects/geowind2")

yrs <- c(1950,2022)

geowind.seas <- read.csv("data/geowind_seasonal_percentiles.csv") %>%
    select(year,season,Q,pc) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")),pc=factor(pc)) %>%
    filter(pc%in%c("5th","50th","95th")) %>%
    na.omit()


nao.eof <- read_table(url("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table"),skip=1,col_names = c("year",1:12))%>%
    gather(month,Value,-year) %>%
    mutate(month=as.numeric(month),
           year=if_else(month==12,year+1,year),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn")) %>%
    group_by(year,season) %>%
    summarise(NAO_EOF=mean(Value))

## nao.gib <- read_table(url("https://crudata.uea.ac.uk/cru/data/nao/nao.dat"),col_names=FALSE,na="-99.99") %>%
##     magrittr::set_colnames(c("year",1:12,"ann")) %>%
##     select(-ann) %>%
##     gather(month,Value,-year) %>%
##     mutate(month=as.numeric(month)) %>%
##         mutate(year=if_else(month==12,year+1,year),
##            season=case_when(month%in%c(12,1,2)~"Winter",
##                             month%in%3:5~"Spring",
##                             month%in%6:8~"Summer",
##                             month%in%9:11~"Autumn")) %>%
##     group_by(year,season) %>%
##     summarise(NAO_Gib=mean(Value))


nao.HUR <- read_table(url("https://climatedataguide.ucar.edu/sites/default/files/2022-10/nao_station_monthly.txt"),col_names=FALSE,na="-999.",skip=2) %>%
    magrittr::set_colnames(c("year",1:12)) %>%
    gather(month,nao,-year) %>%
    mutate(month=as.numeric(month)) %>%
        mutate(month=as.numeric(month)) %>%
        mutate(year=if_else(month==12,year+1,year),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn")) %>%
    group_by(year,season) %>%
    summarise(NAO_HUR=mean(nao))



ea <-read_table("data/ea_index.dat",skip=8) %>%
    rename_all(tolower) %>%
    mutate(year=if_else(month==12,year+1,year),
           season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn")) %>%
    group_by(year,season) %>%
    summarise(EA=mean(index))

cor.table <-geowind.seas %>%
    left_join(nao.HUR) %>%
    left_join(ea) %>%
    left_join(nao.eof) %>%
    filter(between(year,yrs[1],yrs[2]))%>%
    gather(series,tele,-season,-year,-pc,-Q) %>%
    group_by(season,pc,series) %>%
    nest() %>%
    mutate(season=factor(season,levels=c("Winter","Spring","Summer","Autumn")),
           pc=factor(pc,levels=c("95th","50th","5th")),
        model=map(data,~cor.test(.x$Q,.x$tele)),
        tidied=map(model,tidy))%>%
    unnest(tidied) %>%
    select(season,pc,series,estimate,p.value) %>%
    mutate(estimate=round(estimate,2),
           Correlation=if_else(p.value<0.05,paste0(estimate,"*"),as.character(estimate)),
           Correlation=if_else(p.value<0.01,paste0(estimate,"**"),as.character(Correlation))) %>%
    select(Correlation) %>%
    spread(series,Correlation)
    


caption <- sprintf(paste0("Correlation between the geowind percentiles per season and three teleconnection indices over the period %s-%s: EA - East ",
                                "Atlantic Pattern (\\url{https://www.cpc.ncep.noaa.gov/data/teledoc/ea.shtml}); NAO\\_HUR - Hurrell NAO index ",
                                "(\\url{https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based}); ",
                                "NAO\\_EOF - Index Derived from the leading EOF of 500 millibar height field over 0-90N ",
                          "(\\url{https://www.ncdc.noaa.gov/teleconnections/nao/}). * indicates correlations signficant at p$<0.05$; ",
                          " ** indicates correlations significant at p$<0.01$"),yrs[1],yrs[2])

label <- "tab:nao_cor"
print(xtable(cor.table,caption=caption,label=label),file="tables/nao_cor_table.tex",include.rownames=FALSE)
