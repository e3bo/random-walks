#!/usr/bin/env bash

set -e

dvc run \
    -d CEID-InfectionKalman \
    -d CEID-Walk \
    -d COVIDhub-Ensemble \
    -d evaluate-forecasts.R \
    -o figure \
    -o evaluate-forecasts.md \
    -o evaluate-forecasts.html \
    --force \
    -n evaluate \
    'Rscript -e "knitr::spin(\"evaluate-forecasts.R\")"'