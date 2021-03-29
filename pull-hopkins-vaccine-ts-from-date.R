#!/usr/bin/env Rscript

library(magrittr)

fdt <- Sys.getenv("fdt", unset = "2021-03-20")
fdtm <- as.POSIXct(paste(fdt, "0:03:00"), tz = "EST")

if(lubridate::ymd(fdt) < "2021-02-12"){
  print("No vaccine data available for dates before 2021 Feb 12.")
  out <- file.path("hopkins-vaccine", fdt)
  if (!dir.exists(out)){
    dir.create(out)
  }
  q("no")
}

xmlfile <- tempfile()

url <- paste0("https://github.com/govex/COVID-19.git/trunk/data_tables/vaccine_data/us_data/time_series")

cmd <- paste("svn log ", url, " --xml >", xmlfile)
system(cmd)
log <- xml2::read_xml(xmlfile)
datetimes <- 
  xml2::xml_find_all(log, "//logentry/date") %>% 
  xml2::xml_text() %>%
  lubridate::as_datetime()
is_later <- datetimes > fdtm
contemp <- match(FALSE, is_later)
xpath <- paste0("//logentry[", contemp, "]")
revision <- 
  xml2::xml_find_all(log, xpath) %>% 
  xml2::xml_attr("revision")
url2 <- paste0(url, "@r", revision)

if(!dir.exists("hopkins-vaccine")){
  dir.create("hopkins-vaccine")
}
cmd2 <- paste0("svn export ", url2, " ./hopkins-vaccine/", fdt)
system(cmd2)
