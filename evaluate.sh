#!/usr/bin/env bash

set -e

dvc run \
    -d pull-other-forecasts.R \
    -o other-model-forecasts.rds \
    --force \
    -n pull-archived-forecasts \
    ./pull-other-forecasts.R
    
dvc run \
    -d other-model-forecasts.rds \
    -d CEID-InfectionKalman \
    -d make-trajectory-plots.R \
    -o trajectories.png \
    --force \
    -n make-trajectory-plots \
    ./make-trajectory-plots.R