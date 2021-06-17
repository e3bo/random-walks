#!/usr/bin/env Rscript

library(magrittr)

download_past_issue <- function(dname, csvuri) {
  destdir <- file.path("google-mobility-reports-wayback", dname)
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  curdir <- setwd(destdir)
  on.exit(setwd(curdir))
  download.file(url = csvuri, destfile = "Global_Mobility_Report.csv",
                extra = "--retry 3")
  system("head -n1 Global_Mobility_Report.csv > US-states-mobility.csv")
  system("grep \"US-[A-Z][A-Z]\" Global_Mobility_Report.csv >> US-states-mobility.csv")
  unlink("Global_Mobility_Report.csv")
  cat(csvuri, file = "uri")
}

dirname <- Sys.getenv("dname", unset = "2020-06-29")
uri <- Sys.getenv("uri",
                  unset = "https://web.archive.org/web/20200628151643/https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")

download_past_issue(dirname, uri)