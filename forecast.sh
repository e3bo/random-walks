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
    -d pull-healthdata-ts-from-date.R \
    -o healthdata/$fdt/$loc \
    -n healthdata-$fdt-$loc \
    --force \
    fdt=$fdt loc=$loc ./pull-healthdata-ts-from-date.R

dirname=fips-${loc}/${fdt}
mkdir -p $dirname
cd $dirname

dvc run \
    -w ../.. \
    -d InfectionKalmanRegularization.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -o forecasts/${fdt}-fips${loc} \
    --force \
    -n forecast-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript InfectionKalmanRegularization.R