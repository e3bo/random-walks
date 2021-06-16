#!/usr/bin/env Rscript

tstamp <- readr::read_lines("timestamp.txt")
url <- paste0(
  "https://web.archive.org/web/",
  tstamp,
  "/https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv"
)
destpath <- "Global_Mobility_Report.csv"
download.file(url = url, destfile = destpath)