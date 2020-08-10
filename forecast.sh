#!/usr/bin/env bash

set -e

ddt="2020-08-09"
fdt="2020-08-09"
dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$ddt \
    -n hopkins-$ddt \
    --force \
    fdt=$ddt Rscript pull-hopkins-ts-from-date.R
    
dvc run \
    -d forecast-vis.R \
    -d covidhub-common.R \
    -d hopkins/$ddt \
    -o forecasts/$fdt-CEID-Walk.csv \
    --force \
    -n visualize-forecast-$fdt-and-data-$ddt \
    ddt=$ddt fdt=$fdt Rscript forecast-vis.R
