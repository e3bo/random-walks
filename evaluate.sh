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
    -o lambda1000.00-a0.94-CEID-InfectionKalman \
    -o lambda083.33-a0.94-CEID-InfectionKalman \
    -o lambda043.48-a0.94-CEID-InfectionKalman \
    -o lambda029.41-a0.94-CEID-InfectionKalman \
    -o lambda022.22-a0.94-CEID-InfectionKalman \
    -o lambda017.86-a0.94-CEID-InfectionKalman \
    -o lambda014.93-a0.94-CEID-InfectionKalman \
    -o lambda012.82-a0.94-CEID-InfectionKalman \
    -o lambda011.24-a0.94-CEID-InfectionKalman \
    -o lambda010.00-a0.94-CEID-InfectionKalman \
    -o lambda1000.00-a0.95-CEID-InfectionKalman \
    -o lambda083.33-a0.95-CEID-InfectionKalman \
    -o lambda043.48-a0.95-CEID-InfectionKalman \
    -o lambda029.41-a0.95-CEID-InfectionKalman \
    -o lambda022.22-a0.95-CEID-InfectionKalman \
    -o lambda017.86-a0.95-CEID-InfectionKalman \
    -o lambda014.93-a0.95-CEID-InfectionKalman \
    -o lambda012.82-a0.95-CEID-InfectionKalman \
    -o lambda011.24-a0.95-CEID-InfectionKalman \
    -o lambda010.00-a0.95-CEID-InfectionKalman \
    --force \
    ./combine-location-forecasts.R

dvc run \
    -d lambda1000.00-a0.94-CEID-InfectionKalman \
    -d lambda083.33-a0.94-CEID-InfectionKalman \
    -d lambda043.48-a0.94-CEID-InfectionKalman \
    -d lambda029.41-a0.94-CEID-InfectionKalman \
    -d lambda022.22-a0.94-CEID-InfectionKalman \
    -d lambda017.86-a0.94-CEID-InfectionKalman \
    -d lambda014.93-a0.94-CEID-InfectionKalman \
    -d lambda012.82-a0.94-CEID-InfectionKalman \
    -d lambda011.24-a0.94-CEID-InfectionKalman \
    -d lambda010.00-a0.94-CEID-InfectionKalman \
    -d lambda1000.00-a0.95-CEID-InfectionKalman \
    -d lambda083.33-a0.95-CEID-InfectionKalman \
    -d lambda043.48-a0.95-CEID-InfectionKalman \
    -d lambda029.41-a0.95-CEID-InfectionKalman \
    -d lambda022.22-a0.95-CEID-InfectionKalman \
    -d lambda017.86-a0.95-CEID-InfectionKalman \
    -d lambda014.93-a0.95-CEID-InfectionKalman \
    -d lambda012.82-a0.95-CEID-InfectionKalman \
    -d lambda011.24-a0.95-CEID-InfectionKalman \
    -d lambda010.00-a0.95-CEID-InfectionKalman \
    -d calculate-empirical-errors.R \
     -o lambda1000.00-a0.94-CEID-InfectionKalmanEmp \
    -o lambda083.33-a0.94-CEID-InfectionKalmanEmp \
    -o lambda043.48-a0.94-CEID-InfectionKalmanEmp \
    -o lambda029.41-a0.94-CEID-InfectionKalmanEmp \
    -o lambda022.22-a0.94-CEID-InfectionKalmanEmp \
    -o lambda017.86-a0.94-CEID-InfectionKalmanEmp \
    -o lambda014.93-a0.94-CEID-InfectionKalmanEmp \
    -o lambda012.82-a0.94-CEID-InfectionKalmanEmp \
    -o lambda011.24-a0.94-CEID-InfectionKalmanEmp \
    -o lambda010.00-a0.94-CEID-InfectionKalmanEmp \
    -o lambda1000.00-a0.95-CEID-InfectionKalmanEmp \
    -o lambda083.33-a0.95-CEID-InfectionKalmanEmp \
    -o lambda043.48-a0.95-CEID-InfectionKalmanEmp \
    -o lambda029.41-a0.95-CEID-InfectionKalmanEmp \
    -o lambda022.22-a0.95-CEID-InfectionKalmanEmp \
    -o lambda017.86-a0.95-CEID-InfectionKalmanEmp \
    -o lambda014.93-a0.95-CEID-InfectionKalmanEmp \
    -o lambda012.82-a0.95-CEID-InfectionKalmanEmp \
    -o lambda011.24-a0.95-CEID-InfectionKalmanEmp \
    -o lambda010.00-a0.95-CEID-InfectionKalmanEmp \
    --force \
    -n make-empirical-pi-forecasts \
    ./calculate-empirical-errors.R

dvc run \
    -d other-model-forecasts.rds \
    -d lambda1000.00-a0.94-CEID-InfectionKalman \
    -d lambda083.33-a0.94-CEID-InfectionKalman \
    -d lambda043.48-a0.94-CEID-InfectionKalman \
    -d lambda029.41-a0.94-CEID-InfectionKalman \
    -d lambda022.22-a0.94-CEID-InfectionKalman \
    -d lambda017.86-a0.94-CEID-InfectionKalman \
    -d lambda014.93-a0.94-CEID-InfectionKalman \
    -d lambda012.82-a0.94-CEID-InfectionKalman \
    -d lambda011.24-a0.94-CEID-InfectionKalman \
    -d lambda010.00-a0.94-CEID-InfectionKalman \
    -d lambda1000.00-a0.95-CEID-InfectionKalman \
    -d lambda083.33-a0.95-CEID-InfectionKalman \
    -d lambda043.48-a0.95-CEID-InfectionKalman \
    -d lambda029.41-a0.95-CEID-InfectionKalman \
    -d lambda022.22-a0.95-CEID-InfectionKalman \
    -d lambda017.86-a0.95-CEID-InfectionKalman \
    -d lambda014.93-a0.95-CEID-InfectionKalman \
    -d lambda012.82-a0.95-CEID-InfectionKalman \
    -d lambda011.24-a0.95-CEID-InfectionKalman \
    -d lambda010.00-a0.95-CEID-InfectionKalman \
    -d lambda1000.00-a0.94-CEID-InfectionKalmanEmp \
    -d lambda083.33-a0.94-CEID-InfectionKalmanEmp \
    -d lambda043.48-a0.94-CEID-InfectionKalmanEmp \
    -d lambda029.41-a0.94-CEID-InfectionKalmanEmp \
    -d lambda022.22-a0.94-CEID-InfectionKalmanEmp \
    -d lambda017.86-a0.94-CEID-InfectionKalmanEmp \
    -d lambda014.93-a0.94-CEID-InfectionKalmanEmp \
    -d lambda012.82-a0.94-CEID-InfectionKalmanEmp \
    -d lambda011.24-a0.94-CEID-InfectionKalmanEmp \
    -d lambda010.00-a0.94-CEID-InfectionKalmanEmp \
    -d lambda1000.00-a0.95-CEID-InfectionKalmanEmp \
    -d lambda083.33-a0.95-CEID-InfectionKalmanEmp \
    -d lambda043.48-a0.95-CEID-InfectionKalmanEmp \
    -d lambda029.41-a0.95-CEID-InfectionKalmanEmp \
    -d lambda022.22-a0.95-CEID-InfectionKalmanEmp \
    -d lambda017.86-a0.95-CEID-InfectionKalmanEmp \
    -d lambda014.93-a0.95-CEID-InfectionKalmanEmp \
    -d lambda012.82-a0.95-CEID-InfectionKalmanEmp \
    -d lambda011.24-a0.95-CEID-InfectionKalmanEmp \
    -d lambda010.00-a0.95-CEID-InfectionKalmanEmp \
    -d analyze-scores.R \
    -o analyze-scores.md \
    -o analyze-scores.html \
    --plots model.csv \
    --plots location-model.csv \
    --plots horizon-location-model-cases.csv \
    --plots horizon-location-model-hosp.csv \
    --force \
    -n analyze-scores \
    'Rscript -e "knitr::spin(\"analyze-scores.R\")"'

dvc plots modify horizon-location-model-cases.csv --template horizon-location-model
dvc plots modify horizon-location-model-hosp.csv --template horizon-location-model
dvc plots modify location-model.csv --template location-model
dvc plots modify model.csv --template model

dvc run \
    -d other-model-forecasts.rds \
    -d hopkins/$ddt \
    -d healthdata/$ddt \
    -d lambda1000.00-a0.94-CEID-InfectionKalman \
    -d lambda083.33-a0.94-CEID-InfectionKalman \
    -d lambda043.48-a0.94-CEID-InfectionKalman \
    -d lambda029.41-a0.94-CEID-InfectionKalman \
    -d lambda022.22-a0.94-CEID-InfectionKalman \
    -d lambda017.86-a0.94-CEID-InfectionKalman \
    -d lambda014.93-a0.94-CEID-InfectionKalman \
    -d lambda012.82-a0.94-CEID-InfectionKalman \
    -d lambda011.24-a0.94-CEID-InfectionKalman \
    -d lambda010.00-a0.94-CEID-InfectionKalman \
    -d lambda1000.00-a0.95-CEID-InfectionKalman \
    -d lambda083.33-a0.95-CEID-InfectionKalman \
    -d lambda043.48-a0.95-CEID-InfectionKalman \
    -d lambda029.41-a0.95-CEID-InfectionKalman \
    -d lambda022.22-a0.95-CEID-InfectionKalman \
    -d lambda017.86-a0.95-CEID-InfectionKalman \
    -d lambda014.93-a0.95-CEID-InfectionKalman \
    -d lambda012.82-a0.95-CEID-InfectionKalman \
    -d lambda011.24-a0.95-CEID-InfectionKalman \
    -d lambda010.00-a0.95-CEID-InfectionKalman \
    -d make-trajectory-plots.R \
    -o trajectories-all \
    --force \
    -n make-trajectory-plots \
    ./make-trajectory-plots.R