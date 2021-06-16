#!/bin/env bash

dvc run \
  -d timestamp.txt \
  -d download-google-mob.R \
  -o Global_Mobility_Report.csv \
  -n wayback-download \
  --force \
  ./download-google-mob.R