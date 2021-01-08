#!/usr/bin/env bash

set -e

ddt="${ddt:-2020-08-09}"
fdt="${fdt:-2020-08-09}"

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$ddt \
    -n hopkins-$ddt \
    --force \
    fdt=$ddt ./pull-hopkins-ts-from-date.R
    
dvc run \
    -d forecast-vis.R \
    -d covidhub-common.R \
    -d hopkins/$ddt \
    -d forecasts/$fdt-CEID-Walk.csv \
    -O visuals/fdt$fdt-ddt$ddt-inc-death-forecasts.png \
    -O visuals/fdt$fdt-ddt$ddt-cum-death-forecasts.png \
    -O visuals/fdt$fdt-ddt$ddt-inc-case-forecasts.pdf \
    --force \
    -n visualize-forecast-$fdt-and-data-$ddt \
    ddt=$ddt fdt=$fdt ./forecast-vis.R
