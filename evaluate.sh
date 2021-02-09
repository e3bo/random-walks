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
    -o lambda10.8-CEID-InfectionKalman \
    -o lambda15.85-CEID-InfectionKalman \
    -o lambda23.26-CEID-InfectionKalman \
    -o lambda34.15-CEID-InfectionKalman \
    -o lambda50.12-CEID-InfectionKalman \
    -o lambda73.56-CEID-InfectionKalman \
    -o lambda107.98-CEID-InfectionKalman \
    -o lambda158.49-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R
    
dvc run \
    -d other-model-forecasts.rds \
    -d lambda10.8-CEID-InfectionKalman \
    -d lambda15.85-CEID-InfectionKalman \
    -d lambda23.26-CEID-InfectionKalman \
    -d lambda34.15-CEID-InfectionKalman \
    -d lambda50.12-CEID-InfectionKalman \
    -d lambda73.56-CEID-InfectionKalman \
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
    -d lambda10.8-CEID-InfectionKalman \
    -d lambda15.85-CEID-InfectionKalman \
    -d lambda23.26-CEID-InfectionKalman \
    -d lambda34.15-CEID-InfectionKalman \
    -d lambda50.12-CEID-InfectionKalman \
    -d lambda73.56-CEID-InfectionKalman \
    -d lambda107.98-CEID-InfectionKalman \
    -d lambda158.49-CEID-InfectionKalman \
    -d analyze-scores.R \
    -o figure \
    -o analyze-scores.md \
    -o analyze-scores.html \
    --force \
    -n analyze-scores \
    'Rscript -e "knitr::spin(\"analyze-scores.R\")"'