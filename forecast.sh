#!/usr/bin/env bash

set -e

fdt="2020-08-03"
dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$fdt \
    -n hopkins-$fdt \
    --force \
    fdt=$fdt Rscript pull-hopkins-ts-from-date.R
    
dvc run \
    -d forecast-package-covidhub-forecast.R \
    -d model.ini \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -M metrics/$fdt-forecast-calc-time.json \
    -o forecasts/$fdt-CEID-Walk.csv \
    --force \
    -n forecast-$fdt \
    fdt=$fdt Rscript forecast-package-covidhub-forecast.R
