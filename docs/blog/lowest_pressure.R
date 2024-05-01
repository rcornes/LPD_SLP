library(dplyr)
library(dataresqc)
library(lubridate)

## Read the data
data_dir <- "~/projects/geowindWEB/data/SEF"

ifiles <- list.files(data_dir,"London*",full=TRUE)

london <- lapply(ifiles,read_sef,all=TRUE) %>%
    bind_rows()

london %>%
    filter(Value<960) %>%
    mutate(date=make_datetime(Year,Month,Day,Hour,Minute),
           date2=make_date(Year,Month,Day)) %>%
    group_by(date2) %>%
    arrange(Value) %>%
    slice(1) %>%
    ungroup() %>%
    select(date,Value,Meta) %>%
    arrange(Value) %>%
    as.data.frame()
