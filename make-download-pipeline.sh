#!/bin/env bash

set -e
IFS=","
cat wayback-urls.csv | while read dname uri
do
dvc run \
  -d download-google-mob.R \
  -o google-mobility-reports-wayback/$dname \
  -n wayback-download-$dname \
  --force \
  dname=$dname uri=$uri ./download-google-mob.R
done