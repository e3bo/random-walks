#!/usr/bin/env bash

set -e

dvc run \
    -d forecasts \
    -d evaluate-forecasts.R \
    -o figure \
    -o evaluate-forecast.md \
    -o evaluate-forecasts.html \
    --force \
    -n evaluate \
    'Rscript -e "knitr::spin(\"evaluate-forecasts.R\")"'