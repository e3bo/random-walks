#!/usr/bin/env bash

dvc run \
    -d forecast-package-covidhub-forecast.R \
    -d forecast.ini \
    -d model.ini \
    -d covidhub-common.R \
    -d hopkins \
    -M forecast-calc-time.json \
    -o forecasts \
    --force \
    -n forecast \
    Rscript forecast-package-covidhub-forecast.R
