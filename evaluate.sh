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
    -o lambda1000.00-CEID-InfectionKalman \
    -o lambda083.33-CEID-InfectionKalman \
    -o lambda043.48-CEID-InfectionKalman \
    -o lambda029.41-CEID-InfectionKalman \
    -o lambda022.22-CEID-InfectionKalman \
    -o lambda017.86-CEID-InfectionKalman \
    -o lambda014.93-CEID-InfectionKalman \
    -o lambda012.82-CEID-InfectionKalman \
    -o lambda011.24-CEID-InfectionKalman \
    -o lambda010.00-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R

dvc run \
    -d lambda1000.00-CEID-InfectionKalman \
    -d lambda083.33-CEID-InfectionKalman \
    -d lambda043.48-CEID-InfectionKalman \
    -d lambda029.41-CEID-InfectionKalman \
    -d lambda022.22-CEID-InfectionKalman \
    -d lambda017.86-CEID-InfectionKalman \
    -d lambda014.93-CEID-InfectionKalman \
    -d lambda012.82-CEID-InfectionKalman \
    -d lambda011.24-CEID-InfectionKalman \
    -d lambda010.00-CEID-InfectionKalman \
    -d calculate-empirical-errors.R \
    -o lambda1000.00-CEID-InfectionKalmanEmp \
    -o lambda083.33-CEID-InfectionKalmanEmp \
    -o lambda043.48-CEID-InfectionKalmanEmp \
    -o lambda029.41-CEID-InfectionKalmanEmp \
    -o lambda022.22-CEID-InfectionKalmanEmp \
    -o lambda017.86-CEID-InfectionKalmanEmp \
    -o lambda014.93-CEID-InfectionKalmanEmp \
    -o lambda012.82-CEID-InfectionKalmanEmp \
    -o lambda011.24-CEID-InfectionKalmanEmp \
    -o lambda010.00-CEID-InfectionKalmanEmp \
    --force \
    -n make-empirical-pi-forecasts \
    ./calculate-empirical-errors.R

dvc run \
    -d other-model-forecasts.rds \
    -d lambda1000.00-CEID-InfectionKalman \
    -d lambda083.33-CEID-InfectionKalman \
    -d lambda043.48-CEID-InfectionKalman \
    -d lambda029.41-CEID-InfectionKalman \
    -d lambda022.22-CEID-InfectionKalman \
    -d lambda017.86-CEID-InfectionKalman \
    -d lambda014.93-CEID-InfectionKalman \
    -d lambda012.82-CEID-InfectionKalman \
    -d lambda011.24-CEID-InfectionKalman \
    -d lambda010.00-CEID-InfectionKalman \
    -d lambda1000.00-CEID-InfectionKalmanEmp \
    -d lambda083.33-CEID-InfectionKalmanEmp \
    -d lambda043.48-CEID-InfectionKalmanEmp \
    -d lambda029.41-CEID-InfectionKalmanEmp \
    -d lambda022.22-CEID-InfectionKalmanEmp \
    -d lambda017.86-CEID-InfectionKalmanEmp \
    -d lambda014.93-CEID-InfectionKalmanEmp \
    -d lambda012.82-CEID-InfectionKalmanEmp \
    -d lambda011.24-CEID-InfectionKalmanEmp \
    -d lambda010.00-CEID-InfectionKalmanEmp \
    -d analyze-scores.R \
    -o analyze-scores.md \
    -o analyze-scores.html \
    --plots model.csv \
    --plots location-model.csv \
    --plots horizon-location-model.csv \
    --force \
    -n analyze-scores \
    'Rscript -e "knitr::spin(\"analyze-scores.R\")"'

dvc plots modify horizon-location-model.csv --template horizon-location-model    
dvc plots modify location-model.csv --template location-model
dvc plots modify model.csv --template model
    
dvc run \
    -d other-model-forecasts.rds \
    -d lambda1000.00-CEID-InfectionKalmanEmp \
    -d lambda083.33-CEID-InfectionKalmanEmp \
    -d lambda043.48-CEID-InfectionKalmanEmp \
    -d lambda029.41-CEID-InfectionKalmanEmp \
    -d lambda022.22-CEID-InfectionKalmanEmp \
    -d lambda017.86-CEID-InfectionKalmanEmp \
    -d lambda014.93-CEID-InfectionKalmanEmp \
    -d lambda012.82-CEID-InfectionKalmanEmp \
    -d lambda011.24-CEID-InfectionKalmanEmp \
    -d lambda010.00-CEID-InfectionKalmanEmp \
    -d make-trajectory-plots.R \
    -o trajectories-all \
    --force \
    -n make-trajectory-plots \
    ./make-trajectory-plots.R