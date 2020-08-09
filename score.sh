#!/usr/bin/env bash

set -e

ddt="2020-08-03"
dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$ddt \
    -n hopkins-$ddt \
    --force \
    fdt=$ddt Rscript pull-hopkins-ts-from-date.R
    
dvc run \
    -d quantile-score.R \
    -d forecasts \
    -d covidhub-common.R \
    -d hopkins/$ddt \
    --plots-no-cache metrics/$ddt-score-summary.csv \
    -o metrics/$ddt-residuals.rds \
    --force \
    -n score-$ddt \
    ddt=$ddt Rscript quantile-score.R
