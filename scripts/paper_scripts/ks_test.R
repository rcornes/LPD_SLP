library(tidyverse)
library(purrr)
library(broom)
library(xtable)

setwd("~/projects/geowind2")

## Read percentiles
pc_file <- "data/geowind_seasonal_percentiles.csv"

DF <- read.csv(pc_file) %>%
    filter(season=="Winter")

early.rng <- c(1790,1829)
mid.rng <- c(1860,1890)
late.rng <- c(1990,2023)

itest <- map(c("95th","90th","50th","5th"),function(ipc){
    early <- filter(DF,pc==ipc,between(year,early.rng[1],early.rng[2]))
    late <- filter(DF,pc==ipc,between(year,late.rng[1],late.rng[2]))
    mid <- filter(DF,pc==ipc,between(year,mid.rng[1],mid.rng[2]))

    test1 <- ks.test(early$Q,late$Q) %>%
        tidy() %>%
        select("statistic","p.value") %>%
        mutate(pc=ipc,test=sprintf("%s-%s vs. %s-%s",early.rng[1],early.rng[2],late.rng[1],late.rng[2]))

    test2 <- ks.test(mid$Q,late$Q) %>%
        tidy() %>%
        select("statistic","p.value") %>%
        mutate(pc=ipc,test=sprintf("%s-%s vs. %s-%s",mid.rng[1],mid.rng[2],late.rng[1],late.rng[2]))

    test3 <- ks.test(early$Q,mid$Q) %>%
        tidy() %>%
        select("statistic","p.value") %>%
        mutate(pc=ipc,test=sprintf("%s-%s vs. %s-%s",early.rng[1],early.rng[2],mid.rng[1],mid.rng[2]))

    op <- bind_rows(test1,test2,test3) %>%
        mutate(p.value=round(p.value,4)) %>%
        rename(Percentile=pc,`Test Periods`=test)

    op
})

itest <- bind_rows(itest)
caption=paste("Results from Kolmogorov-Smirnov tests (two-sample, two-sided) that evaluate geowind percentiles over different periods.",
               "Note that the results are sensitive to the exact years that are contained within each sample.",
               "In certain comparisons ties are present in the data, which preclude the calculation of exact p-values.")
print(xtable(itest,caption=caption,label="tab:kstest"),file="tables/KS_test.tex",include.rownames=FALSE)
