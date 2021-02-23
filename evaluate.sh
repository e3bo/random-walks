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
    -d analyze-scores.R \
    -o figure \
    -o analyze-scores.md \
    -o analyze-scores.html \
    --plots wis-model-location.csv \
    --plots wis-horizon-location-model.csv \
    --force \
    -n analyze-scores \
    'Rscript -e "knitr::spin(\"analyze-scores.R\")"'

dvc plots modify wis-horizon-location-model.csv --template relscore-horizon-location-model    
dvc plots modify wis-model-location.csv --template relscore-location-model