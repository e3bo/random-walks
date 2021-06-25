#!/usr/bin/env bash

set -e

ddt="${ddt:-2021-02-01}"

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
    -o lambda020.00-status-quo-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R

dvc run \
    -d lambda020.00-status-quo-CEID-InfectionKalman \
    -d other-model-forecasts.rds \
    -d analyze-scores.R \
    -o hosp-wis-by-date.png \
    -o cases-deaths-wis-by-date.png \
    -o cases-deaths-wis.png \
    -o hosp-wis.png \
    --force \
    -n analyze-scores \
    ./analyze-scores.R

dvc run \
    -d other-model-forecasts.rds \
    -d hopkins/$ddt \
    -d healthdata/$ddt \
    -d lambda020.00-status-quo-CEID-InfectionKalman \
    -d make-trajectory-plots.R \
    -o trajectories-all \
    --force \
    -n make-trajectory-plots \
    ./make-trajectory-plots.R