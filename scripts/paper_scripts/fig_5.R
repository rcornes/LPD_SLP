#!/usr/bin/env Rscript

## --------------------------------------------------------------------------------
## Script to product the plot for Figure 4 and correlation table of
## jet parameters and geostrophic wind quantiles
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-04
##
## Pre-requisites: download_20CR_jet.sh, make_percentiles.R
## --------------------------------------------------------------------------------

correlation <- function(x,y,detrend=TRUE,...){
    require(broom)
    require(dplyr)
    
    C <- cor.test(x,y,...)
    stats <- tidy(C) %>%
        mutate(series="correlation")

    
    if(detrend){
        x.detrend <- residuals(lm(x~c(1:length(x))))
        y.detrend <- residuals(lm(y~c(1:length(y))))
        C2 <- cor.test(x.detrend,y.detrend,...)
        stats <- tidy(C2) %>%
            mutate(series="correlation detrended") %>%
            rbind(stats)
    }

    stats
    
}

lanczos_weights<-function(window=101,sampl_rate=1,type="lowpass",
                          low_freq=1/100,high_freq=1/10){
    
    ## Code taken from: https://stackoverflow.com/questions/17264119/using-lanczos-low-pass-filter-in-r-program
    ## Based on the notes available here: https://www2.atmos.umd.edu/~ekalnay/syllabi/AOSC630/METO630ClassNotes13.pdf
    
    low_freq <- sampl_rate*low_freq
    high_freq <- sampl_rate*high_freq

    if (type=="lowpass"){
        order <- ((window - 1) %/% 2 ) + 1
        nwts <- 2 * order + 1
        fc <- low_freq
        w <- seq(0,0,length=nwts)
        n <- nwts %/% 2
        w[n+1] <- 2 * fc
        k <- seq(1, n-1)
        sigma <- sin(pi * k / n) * n / (pi * k)
        firstfactor <- sin(2 *pi * fc * k) / (pi * k)
        w[n:2] <- firstfactor * sigma
        w[(n+2):(length(w)-1)] <- firstfactor * sigma
        w <- w[-c(1,length(w))]}
    else if (type=="highpass"){
            order <- ((window - 1) %/% 2 ) + 1
            nwts <-  2 * order + 1
            fc <- high_freq
            w  <- seq(0,0,length=nwts)
            n  <- nwts %/% 2
            w[n+1] <- 2 * fc
            k <- seq(1, n-1)
            sigma <- sin(pi * k / n) * n / (pi * k)
            firstfactor <- sin(2 *pi * fc * k) / (pi * k)
            w[n:2] <- firstfactor * sigma
            w[(n+2):(length(w)-1)] <- firstfactor * sigma
            w <- w[-c(1,length(w))]
            w <- -w
            w[order] <- 1-2*fc }
    else if (type=="bandpass"){
                order <- ((window - 1) %/% 2 ) + 1
                nwts <- 2 * order + 1
                fc <- low_freq
                w <- seq(0,0,length=nwts)
                n <- nwts %/% 2
                w[n+1] <- 2 * fc
                k <- seq(1, n-1)
                sigma <- sin(pi * k / n) * n / (pi * k)
                firstfactor <- sin(2 *pi * fc * k) / (pi * k)
                w[n:2] <- firstfactor * sigma
                w[(n+2):(length(w)-1)] <- firstfactor * sigma
                w1 <- w[-c(1,length(w))]

                order <- ((window - 1) %/% 2 ) + 1
                nwts <- 2 * order + 1
                fc <- high_freq
                w <- seq(0,0,length=nwts)
                n <- nwts %/% 2
                w[n+1] <- 2 * fc
                k <- seq(1, n-1)
                sigma <- sin(pi * k / n) * n / (pi * k)
                firstfactor <- sin(2 *pi * fc * k) / (pi * k)
                w[n:2] <- firstfactor * sigma
                w[(n+2):(length(w)-1)] <- firstfactor * sigma
                w2 <- w[-c(1,length(w))]
                w <- w2-w1}
    else {print("Please specify a valid filter type: either 'lowpass', 'highpass' or 'bandpass'")}
    return(w)
}

library(plyr)
library(tidyverse)
library(lubridate)
library(broom)
library(ncdf4)
library(rnaturalearth)
library(rnaturalearthdata)
library(lattice)
library(latticeExtra)
library(magrittr)
library(ggthemes)
library(zyp)
library(gridExtra)
library(cowplot)
library(xtable)
library(metR)
library(tactile)

setwd("~/projects/geowind2")

rng1 <- c(1871,1989)
rng2 <- c(1990,2015)

pc_file <- "data/geowind_seasonal_percentiles.csv"
geodat <- read.csv(pc_file) %>%
    filter(pc%in%c("5th","50th","95th")) %>%
    mutate(season=factor(season,levels=c("Summer","Autumn","Winter","Spring")),
           pc=factor(pc)) %>%
    na.omit()

## Calculate jet parameters
jet.file <- nc_open("data/20CRv3_zonmean.nc")
uwnd <- ncvar_get(jet.file)
uwnd.time <- hours(ncvar_get(jet.file,"time"))+make_datetime(1800)
uwnd.lat <- ncvar_get(jet.file,"lat")
nc_close(jet.file)

filter.wgts <- lanczos_weights(window=61,low_freq=1/10)
class(filter.wgts) <- "Ma"

## Jet Stats - Ensemble mean
uwnd.filt <- apply(uwnd,1,function(x) signal::filtfilt(filter.wgts,x))
jet.speed <- apply(uwnd.filt,1,max)
jet.lat <- apply(uwnd.filt,1,function(x) uwnd.lat[which.max(x)])

jet.seas <- data.frame(time=uwnd.time,speed=jet.speed,lat=jet.lat) %>%
    gather(series,data,-time) %>%
        mutate(month=month(time),year=year(time),
                      season=case_when(month%in%c(12,1,2)~"Winter",
                            month%in%3:5~"Spring",
                            month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")),
           year=ifelse(month%in%12,year+1,year)) %>%
    group_by(year,season,series) %>%
    summarise(av=mean(data)) %>%
    spread(series,av)

## Jet Stats - Ensemble members
jet.ens.file <- nc_open("data/20CRv3_ens_zonmean.nc")
uwnd.ens.time <- hours(floor(ncvar_get(jet.ens.file,"time")))+make_datetime(1800)
uwnd.ens.lat <- ncvar_get(jet.ens.file,"lat")

nens <- 80

ens.list <- lapply(1:nens,function(i){
    print(i)
    uwnd <- ncvar_get(jet.ens.file,start=c(1,1,i,1),count = c(-1,-1,1,-1))
    uwnd.filt <- apply(uwnd,1,function(x) signal::filtfilt(filter.wgts,x))
    jet.speed <- apply(uwnd.filt,1,max)
    jet.lat <- apply(uwnd.filt,1,function(x) uwnd.ens.lat[which.max(x)])

    jet.seas <- data.frame(time=uwnd.ens.time,speed=jet.speed,lat=jet.lat) %>%
        tibble() %>%
        gather(series,data,-time) %>%
        mutate(month=month(time),year=year(time),
               season=case_when(month%in%c(12,1,2)~"Winter",
                                month%in%3:5~"Spring",
                                month%in%6:8~"Summer",
                            month%in%9:11~"Autumn"),
           season=factor(season,levels=c("Winter","Spring","Summer","Autumn")),
           year=ifelse(month%in%12,year+1,year)) %>%
    group_by(year,season,series) %>%
    summarise(av=mean(data))

    jet.seas
})

nc_close(jet.ens.file)

jet.ens.seas <- do.call("rbind",ens.list) %>%
    group_by(year,season,series) %>%
    summarise(sd=sd(av))

df.djf <- geodat %>%
    inner_join(jet.seas) %>%
    filter(between(year,rng1[1],rng2[2]),season=="Winter") %>%
    select(year,Q,speed,lat,pc)

## Correlation table between geowind and jet parameters (Supp. Info.)
cor.table <- df.djf %>%
    gather(jet_stat,jet,-year,-Q,-pc) %>%
    group_by(jet_stat,pc) %>%
    nest() %>%
    mutate(cor=map(data,~correlation(.x$Q,.x$jet))) %>%
    unnest(cor) %>%
    select(estimate,p.value,series) %>%
    mutate_at(vars(estimate),round,digits=2) %>%
    mutate_at(vars(p.value),round,digits=4)

caption <- "Correlations between Jet Parameters and Geowind Statistics"
label <- "tab:jetcor"
print(xtable(cor.table,caption=caption,label=label),file="tables/jet_correlations.tex",include.rownames=FALSE)


## Jet/geowind time series plots
plot.df <- df.djf %>%
    select(year,Q,speed,lat,pc) %>%
    spread(pc,Q) %>%
    gather(series,Q,-year) %>%
    filter(series%in%c("lat","speed","95th","50th")) %>%
    left_join(jet.ens.seas) %>%
    mutate(lo=Q-(2*sd),hi=Q+(2*sd),lo=ifelse(is.na(lo),0,lo),hi=ifelse(is.na(hi),0,hi)) %>%
    mutate(series=factor(series,levels=c("50th","95th","speed","lat"),
                         labels=c("Geowind (50th)","Geowind (95th)",
                                  "Jet Speed","Jet Latitude")))

cols <- c("#0080ff","#ff00ff","darkgreen","#ff0000")


plot.df.lat <- filter(plot.df,series=="Jet Latitude")
ts.plot.lat <- xyplot(Q~year|series,data=plot.df.lat,
                       group=series,type="l",
                       scales=list(tck=c(1,0,1,1),y=list(relation="free",cex=0.5),
                                   x=list(cex=0.5)),par.strip.text=list(cex=0.5),
                       xlab=NULL,ylab=list(label=substitute(paste("Latitude (",degree,"N)")),
                                           cex=0.5), lwd=0.75,
                       par.settings=list(layout.heights=list(main.key.padding=0,
                                                             key.axis.padding=0,
                                                             bottom.padding=-1.25,
                                                             top.padding=1),
                                         par.main.text = list(cex=0.75,font = 0.5,
                                                              just = "left",
                                                              x = grid::unit(5, "mm"))),
                       xscale.components = xscale.components.subticks,
                      yscale.components = yscale.components.subticks,
                      lower=plot.df.lat$lo,upper=plot.df.lat$hi,
                       panel=function(x, y,...){
                           col=cols[4]
                           M <- zyp.zhang(y,x)
                           trend <- sprintf("%3.2f",round(M[["trend"]] * 100,2))
                           ## lo.trend <- sprintf("%3.2f",round(M[["lbound"]]*10,2))
                           ## hi.trend <- sprintf("%3.2f",round(M[["ubound"]]*10,2))
                           star <- if(M[["sig"]]<0.05) "**" else if(M[["sig"]]<0.1) "*" else ""
                           if(panel.number()==4) {
                               label <- substitute(paste("Trend = ",trend,star,
                                                         degree,"lat ", century^-1))
                      } else {
                          label <- substitute(paste("Trend = ",trend,star,degree," ", century^-1))
                      }
                           panel.ci(x,y,alpha=0.4,...)
                      panel.xyplot(x,y,col=col,...)
                      panel.abline(b=M[["trend"]],a=M[["intercept"]],col=col)
                      panel.text(1865,max(y)-max(y)*0.01,col=col,labels=label,pos=4,cex=0.5) })

plot.df.speed <- filter(plot.df,!series=="Jet Latitude")

ts.plot.speed <- xyplot(Q~year|series,data=plot.df.speed,
                         group=series,type="l",
                         layout=c(1,3),scales=list(tck=c(1,0,1,1),y=list(relation="free",cex=0.5),
                                                   x=list(cex=0.5)),par.strip.text=list(cex=0.5),
                         xlab=list(label="Year",cex=0.5),ylab=list(label=substitute(paste("Geostrophic Wind Speed (m ",s^-1,")")),
                                                                   cex=0.5),lwd=0.75,
                         par.settings=list(layout.heights=list(top.padding=-5),
                                           par.main.text = list(cex=0.75,font = 0.5,
                                                                just = "left",
                                                                x = grid::unit(5, "mm"))),
                         xscale.components = xscale.components.subticks,
                        yscale.components = yscale.components.subticks,
                        lower=plot.df.speed$lo,
                        upper=plot.df.speed$hi,
                         panel=function(x, y,...){
                             col=cols[panel.number()]
                             M <- zyp.zhang(y,x)
                             trend <- sprintf("%3.2f",round(M[["trend"]] * 100,2))
                             #lo.trend <- sprintf("%3.2f",round(M[["lbound"]]*10,2))
                                        #hi.trend <- sprintf("%3.2f",round(M[["ubound"]]*10,2))
                             star <- if(M[["sig"]]<0.05) "**" else if(M[["sig"]]<0.1) "*" else ""
                             if(panel.number()==4) {
                                 label <- substitute(paste("Trend = ",trend,star,
                                                           degree,"lat ", century^-1))
                             } else {
                                 label <- substitute(paste("Trend = ",trend,star," m ",s^-1," ",
                                                           century^-1))
                             }
                             panel.ci(x,y,alpha=0.4,...)
                             panel.xyplot(x,y,col=col,...)
                             panel.abline(b=M[["trend"]],a=M[["intercept"]],col=col)
                             panel.text(1865,max(y)-max(y)*0.01,col=col,labels=label,
                                        pos=4,cex=0.5) })

ts.plot <- arrangeGrob(ts.plot.lat,ts.plot.speed,heights=c(1.65/4,3/4))

## Read 20CR U850 data
nc.u850 <- nc_open("data/20CR_v3_U850_NA_seas.nc")
u850 <- ncvar_get(nc.u850)
u850.lat <- ncvar_get(nc.u850,"lat")
u850.lon <- ncvar_get(nc.u850,"lon")
u850.time <- ncvar_get(nc.u850,"time")
u850.time <- hours(u850.time)+make_datetime(1800)
dimnames(u850) <- list(u850.lon,u850.lat,u850.time)
nc_close(nc.u850)

## Averages of jet params over date ranges
ind1 <- which(month(u850.time)==1&year(u850.time)>=rng1[1]&year(u850.time)<=rng1[2])
ind2 <- which(month(u850.time)==1&year(u850.time)>=rng2[1]&year(u850.time)<=rng2[2])

u850.early <- apply(u850[,,ind1],c(1,2),FUN=mean)
u850.late <- apply(u850[,,ind2],c(1,2),FUN=mean)

u850.early <- adply(u850.early,c(1,2)) %>%
    as_tibble() %>%
    set_colnames(c("lon","lat","data"))%>%
    mutate(lat=as.numeric(as.character(lat)),lon=as.numeric(as.character(lon))) %>%
    mutate(lon=if_else(lon>180,-360+lon,lon)) %>%
    filter(between(lat,15,85))

u850.late <- adply(u850.late,c(1,2)) %>%
    as_tibble() %>%
    set_colnames(c("lon","lat","data"))%>%
    mutate(lat=as.numeric(as.character(lat)),lon=as.numeric(as.character(lon))) %>%
    mutate(lon=if_else(lon>180,-360+lon,lon)) %>%
    filter(between(lat,15,85))

u850.line.df <- bind_rows(mutate(u850.early,period=sprintf("%s-%s",rng1[1],rng1[2])),
                          mutate(u850.late,period=sprintf("%s-%s",rng2[1],rng2[2])))

## Jet Map Plot
world <- ne_coastline(scale = "medium", returnclass = "sf")
lon.lim <- c(-60,20)
lat.lim <- c(25,70)

u850.map <- ggplot()+
     geom_contour_fill(aes(lon,lat,z=data,fill=after_stat(level)),
                       data=filter(u850.line.df,data>=4),breaks = seq(6,14,1))+
     geom_contour(aes(lon,lat,z=data),data=filter(u850.line.df,data>=4),
                  breaks = seq(6,14,1),lwd=0.1,color="black")+
     scale_fill_brewer(palette = "YlOrRd")+
     geom_sf(data=world,size=0.25,lwd=0.08)+
     geom_text_contour(aes(lon,lat,z=data),data=filter(u850.line.df,data>=4),
                       stroke = 0.2,breaks=seq(6,14,1),size=1,skip=0)+
     theme_map(8)+
     labs(fill=expression(paste("U850 (m/s)")),linetype="U850 average")+
     guides(linetype=guide_legend(title.position="top",direction="vertical"),
            fill=guide_colourbar(title.position="top"))+
     xlim(lon.lim)+
     ylim(lat.lim)+
     theme(plot.title=element_text(face="plain"),legend.position = "bottom",
           legend.key.height=unit(0.2,"cm"),legend.key.width=unit(0.5,"cm"),
           plot.background=element_rect(fill="white", colour="white"),
           rect = element_rect(fill = 'grey'),
           legend.text=element_text(size=4),
           legend.title=element_text(size=4),
           strip.text = element_text(size = 5.5),
           plot.margin = margin(l=30,t=15,b=1,r=1),
           strip.background = element_rect(fill="#ffe5cc"))+
     coord_sf(xlim=c(lon.lim[1]+4,lon.lim[2]-4),ylim=c(lat.lim[1]+4,lat.lim[2]-4))+
     facet_wrap(~period,ncol=2)

pdf("plots/fig_5.pdf",width=3.5,height=6)
plot_grid(u850.map,ts.plot,ncol=1,rel_heights=c(0.56,2),rel_widths = c(0.56,1),
          labels = c("a)","b)"),label_size = 8)
dev.off()
