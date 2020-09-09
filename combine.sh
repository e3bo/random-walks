#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-08-03}"

dvc run \
    -d drift/forecasts/$fdt-CEID-Walk.csv \
    -d no-drift/forecasts/$fdt-CEID-Walk.csv \
    -d combine-forecasts.R \
    -d model-weights/$fdt-weights.csv \
    --force \
    -n combine-$fdt \
    -o forecasts/$fdt-CEID-Walk.csv \
    fdt=$fdt ./combine-forecasts.R
