#!/usr/bin/env bash

set -e

dvc run \
    -d pull-other-forecasts.R \
    -o other-model-forecasts.rds \
    --force \
    -n pull-archived-forecasts \
    ./pull-other-forecasts.R

dvc run \
    -d forecasts \
    -d covidhub-common.R \
    -d combine-location-forecasts.R \
    -n combine \
    -o lambda010.80-CEID-InfectionKalman \
    -o lambda015.85-CEID-InfectionKalman \
    -o lambda023.26-CEID-InfectionKalman \
    -o lambda034.15-CEID-InfectionKalman \
    -o lambda050.12-CEID-InfectionKalman \
    -o lambda073.56-CEID-InfectionKalman \
    -o lambda107.98-CEID-InfectionKalman \
    -o lambda158.49-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R
    
dvc run \
    -d other-model-forecasts.rds \
    -d lambda010.80-CEID-InfectionKalman \
    -d lambda015.85-CEID-InfectionKalman \
    -d lambda023.26-CEID-InfectionKalman \
    -d lambda034.15-CEID-InfectionKalman \
    -d lambda050.12-CEID-InfectionKalman \
    -d lambda073.56-CEID-InfectionKalman \
    -d lambda107.98-CEID-InfectionKalman \
    -d lambda158.49-CEID-InfectionKalman \
    -d make-trajectory-plots.R \
    -o trajectories-all.png \
    -o trajectories-0.png \
    -o trajectories-1.png \
    -o trajectories-2.png \
    -o trajectories-3.png \
    --force \
    -n make-trajectory-plots \
    ./make-trajectory-plots.R
    
dvc run \
    -d other-model-forecasts.rds \
    -d lambda010.80-CEID-InfectionKalman \
    -d lambda015.85-CEID-InfectionKalman \
    -d lambda023.26-CEID-InfectionKalman \
    -d lambda034.15-CEID-InfectionKalman \
    -d lambda050.12-CEID-InfectionKalman \
    -d lambda073.56-CEID-InfectionKalman \
    -d lambda107.98-CEID-InfectionKalman \
    -d lambda158.49-CEID-InfectionKalman \
    -d analyze-scores.R \
    -o figure \
    -o analyze-scores.md \
    -o analyze-scores.html \
    --force \
    -n analyze-scores \
    'Rscript -e "knitr::spin(\"analyze-scores.R\")"'