#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-07-26}"
loc="36"

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$fdt \
    -n hopkins-$fdt \
    --force \
    fdt=$fdt ./pull-hopkins-ts-from-date.R
    
dvc run \
    -d infection-kalman-data-prep.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -o data--$fdt--$loc.csv \
    -o initial-pars--$fdt--$loc.csv \
    --force \
    -n data-prep-$fdt-$loc \
    fdt=$fdt Rscript infection-kalman-data-prep.R
