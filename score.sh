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
    -d forecasts \
    -d covidhub-common.R \
    -d combine-location-forecasts.R \
    -n combine \
    -o lambda125.89-CEID-InfectionKalman \
    -o lambda158.49-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R