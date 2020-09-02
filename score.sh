#!/usr/bin/env bash

set -e

ddt="${ddt:-2020-08-03}"

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$ddt \
    -n hopkins-$ddt \
    --force \
    fdt=$ddt ./pull-hopkins-ts-from-date.R
    
dvc run \
    -d quantile-score.R \
    -d forecasts \
    -d covidhub-common.R \
    -d hopkins/$ddt \
    --plots-no-cache metrics/$ddt-score-by-loc-type.csv \
    --plots-no-cache metrics/$ddt-score-by-loc-type-targ-type.csv \
    --plots-no-cache metrics/$ddt-score-by-loc-type-targ-type-forecast-date.csv \
    -o metrics/$ddt-residuals.rds \
    --force \
    -n score-$ddt \
    ddt=$ddt ./quantile-score.R
dvc plots modify -t lqs -x loc_type --x-label "Location type" metrics/${ddt}-score-by-loc-type.csv
dvc plots modify -t lqs-loctype-panel -x target_type metrics/${ddt}-score-by-loc-type-targ-type.csv
dvc plots modify -t lqs-loctype-targtype-panels -x forecast_date metrics/${ddt}-score-by-loc-type-targ-type-forecast-date.csv
