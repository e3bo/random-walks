#!/usr/bin/env Rscript

library(tidyverse)
library(covidcast)

fdt <- Sys.getenv("fdt", unset = "2021-03-29")
out <- file.path("covidcast-safegraph-home-prop-7dav", fdt)
if(!dir.exists(out)){
  dir.create(out, recursive = TRUE)
}

if(lubridate::ymd(fdt) < "2019-01-01"){
  print("No safegraph data available for dates before 2019 Jan 01.")
  q("no")
}

df <- covidcast_signal(data_source = "safegraph", 
                 signal = "completely_home_prop_7dav",
                 start_day = "2019-01-01", end_day = fdt,
                 geo_type = "state", as_of = fdt)

fp <- file.path(out, "epidata.csv")
write_csv(df, path = fp)