#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-07-26}"
loc="${loc:-36}"

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$fdt \
    -n hopkins-$fdt \
    --force \
    fdt=$fdt ./pull-hopkins-ts-from-date.R

dvc run \
    -d InfectionKalmanRegularization.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -o forecasts/${fdt}-fips${loc}-lambda158.49-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda107.98-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda073.56-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda050.12-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda034.15-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda023.26-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda015.85-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda010.80-CEID-InfectionKalman.csv \
    --force \
    -n forecast-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript InfectionKalmanRegularization.R