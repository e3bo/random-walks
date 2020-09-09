#!/usr/bin/env bash

set -e

ddt="${ddt:-2020-08-03}"

dvc run \
    -d drift/metrics/$ddt-residuals.rds \
    -d no-drift/metrics/$ddt-residuals.rds \
    -d optimally-weight.R \
    -o model-weights/$ddt-weights.csv \
    --force \
    -n weight-$ddt \
    ddt=$ddt ./optimally-weight.R
