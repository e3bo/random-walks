#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-07-26}"

dvc run \
    -d ../forecast-package-covidhub-forecast.R \
    -d model.ini \
    -d ../covidhub-common.R \
    -d ../hopkins/$fdt \
    -M metrics/$fdt-forecast-calc-time.json \
    -o forecasts/$fdt-CEID-Walk.csv \
    --force \
    -n forecast-$fdt \
    fdt=$fdt Rscript ../forecast-package-covidhub-forecast.R
