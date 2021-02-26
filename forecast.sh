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

dirname=fips-${loc}
mkdir -p $dirname
cd $dirname

dvc run \
    -w .. \
    -d InfectionKalmanRegularization.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -o forecasts/${fdt}-fips${loc}-lambda1000.00-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda083.33-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda043.48-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda029.41-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda022.22-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda017.86-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda014.93-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda012.82-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda011.24-CEID-InfectionKalman.csv \
    -o forecasts/${fdt}-fips${loc}-lambda010.00-CEID-InfectionKalman.csv \
    --force \
    -n forecast-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript InfectionKalmanRegularization.R