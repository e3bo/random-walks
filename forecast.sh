#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-07-26}"
loc="36"

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$fdt \
    -n hopkins-$fdt \
    --force \
    fdt=$fdt ./pull-hopkins-ts-from-date.R
    
dvc run \
    -d infection-kalman-data-prep.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -o data--$fdt--$loc.csv \
    -o initial-pars--$fdt--$loc.csv \
    --force \
    -n data-prep-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript infection-kalman-data-prep.R
    
dvc run \
    -d data--$fdt--$loc.csv \
    -d initial-pars--$fdt--$loc.csv \
    -d InfectionKalman.jl \
    -d InfectionKalmanMain.jl \
    -o minimizer--$fdt--$loc.csv \
    --force \
    -n data-fit-$fdt-$loc \
    julia InfectionKalmanMain.jl $fdt $loc
    
dvc run \
    -d minimizer--$fdt--$loc.csv \
    -d data--$fdt--$loc.csv \
    -d infection-kalman-process-results.R \
    -o forecasts/$fdt-CEID-InfectionKalman.csv \
    --force \
    -n write-forecast-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript infection-kalman-process-results.R
