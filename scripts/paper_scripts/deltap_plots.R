#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to produce the Supp. Info. station-based storm indices
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Prerequisites scripts: make_geowind.R
## --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(broom)

setwd("~/projects/geowind2")

source("scripts/filter_CRU.r")
yrs <- c(1748,2023)

geodat <- readRDS("data/geodat_QC.RDS") %>%
    filter(minute(time)==0,hour(time)==12,between(year(time),yrs[1],yrs[2])) %>%
    select(time,debilt,london,paris, contains("QC")) %>%
    mutate(london=if_else(!geowindQC==0,NA_real_,london),
           paris=if_else(!geowindQC==0,NA_real_,paris),
           debilt=if_else(!geowindQC==0,NA_real_,debilt),
           london.lag=london-lag(london), paris.lag=paris-lag(paris), debilt.lag=debilt-lag(debilt))

geodat.seas <- geodat %>%
    mutate(month=month(time),year=year(time),
           season=case_when(month%in%c(12,1,2)~"DJF",
                            month%in%3:5~"MAM",
                            month%in%6:8~"JJA",
                            month%in%9:11~"SON"),
           year=ifelse(month%in%12,year+1,year)) %>%
    select(year,season,london,paris,debilt,contains("lag")) %>%
    gather(series,data,-year,-season)

geodat.lseas <- geodat %>%
    mutate(month=month(time),year=year(time),
           season=case_when(month%in%c(12,1,2,3)~"DJFM"),
           year=ifelse(month%in%12,year+1,year)) %>%
    filter(season=="DJFM") %>%
    select(year,season,london,paris,debilt,contains("lag")) %>%
    gather(series,data,-year,-season)

geodat.ann <- geodat %>%
    mutate(month=month(time), year=year(time), season="Annual") %>%
    select(year,season,london,paris,debilt,contains("lag")) %>%
    gather(series,data,-year,-season)

df <- rbind(geodat.seas,geodat.lseas,geodat.ann)


press.stats <- function(x,thresh=0.3){
    stats <- tibble(n1005=length(which(x<1005)), n980=length(which(x<980)),q1=quantile(x,0.1,na.rm=TRUE))

    miss <- length(which(is.na(x)))/length(x)
    if(miss>thresh) stats[] <- NA_real_

    stats
}

press.diff.stats <- function(x,thresh=0.3){
    stats <- tibble(absdiff_16=length(which(abs(x)>16)), p95diff=quantile(x,0.95,na.rm=TRUE), p99diff=quantile(x,0.99,na.rm=TRUE),
                    mndiff=mean(abs(x),na.rm=TRUE))

    miss <- length(which(is.na(x)))/length(x)
    if(miss>thresh) stats[] <- NA_real_

    stats
}

deltap <- df %>%
    filter(!grepl("lag",series)) %>%
    group_by(year,season,series) %>%
    summarise(press.stats(data)) %>%
    gather(stat,data,-year:-series)

deltap.diff <- df %>%
    filter(grepl("lag",series))%>%
    group_by(year,season,series) %>%
    summarise(press.diff.stats(data)) %>%
    gather(stat,data,-year:-series) %>%
    separate(series,c("series",NA))

istats <- rbind(deltap,deltap.diff) %>%
    mutate(season=factor(season,levels=c("Annual","DJFM","DJF","MAM","JJA","SON")))

## PCA
istats_pca <- istats %>%
    spread(stat,data) %>%
    select(-n980) %>%
    na.omit() %>%
    group_by(series,season) %>%
    nest(data = c(year, absdiff_16, mndiff, n1005, p95diff, p99diff, q1)) %>%
    mutate(pca = map(data, ~ prcomp(.x %>% dplyr::select(-year), center = TRUE, scale = TRUE)),
           pca_aug = map2(pca, data, ~augment(.x, data = .y)))

istats_pca_loadings <- istats_pca %>%
    select(season,series,pca) %>%
    group_by(season,series) %>%
    mutate(loadings=map(pca,tidy,matrix="loadings")) %>%
    unnest(loadings) %>%
    filter(PC==1) %>%
    select(-pca) %>%
    as.data.frame()

## Calculate the polarity for PC1 (assuming these are the same across the different indices)
PC1_multiplier <- istats_pca_loadings %>%
    filter(column=="absdiff_16") %>%
    mutate(mult=if_else(value<0,-1,1)) %>%
    select(season,series,mult)

var_exp <- istats_pca %>% 
    unnest(pca_aug) %>%
    group_by(series,season) %>%
    summarize_at(.vars = vars(contains(".fittedPC")), .funs = funs(var)) %>% 
    gather(key = pc, value = variance,-season,-series) %>%
    group_by(series,season) %>%
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", "")) %>%
  as.data.frame()

## The multuplier ensures that positive values in loadings == more stormy
PC1 <- istats_pca %>%
    unnest(pca_aug) %>%
    select(year,season,series,.fittedPC1) %>%
    rename(data=.fittedPC1) %>%
    mutate(stat="PC1") %>%
    left_join(PC1_multiplier) %>%
    mutate(data=data*mult) %>%
    select(-mult)

istats <- rbind(istats,PC1)

## Gaussian filter
istats <- istats %>%
    group_by(season,series,stat) %>%
    mutate(smooth=filter.cru(thalf=51,tsin=data)[[1]])

storm.stats <- unique(istats$stat)
storm.labs <- paste0(letters[1:length(storm.stats)],") ", storm.stats)
ylabs <- list("n1005"=bquote(N[1005]), "n980"=bquote(N[980]),"q1"=bquote(P[10] (hPa)),"absdiff_16"=bquote(N[Delta~p/Delta~t]),
              "p95diff"=bquote(P[95]*Delta~p),"p99diff"=bquote(P[99]*Delta~p),"mndiff"=bquote(bar(Delta~p/Delta~t)),"PC1"="PC1")

iplots <- map(seq_along(storm.labs),function(i){
    istat <- storm.stats[i]
    
    G <- istats %>%
    filter(between(year,1748,2023)) %>%
    filter(stat==istat) %>%
    ggplot(aes(year,data))+
        facet_grid(season~series,scales="free_y")+
        labs(title=storm.labs[i],x="Year",y=ylabs[[istat]])+
        theme_bw()

    if(istat%in%c("n1005","n980","absdiff_16")){
        G <- G+geom_step()
    } else {
        G <- G + geom_line()
    }

    if(istat=="q1") G <- G + scale_y_reverse()

    G+geom_line(aes(y=smooth),col="red")

})

isave <- map2(storm.stats, iplots, function(L,P) ggsave(file.path("plots",sprintf("deltap_%s.pdf",L)),P))
    
