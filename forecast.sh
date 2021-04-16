#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-12-07}"
loc="${loc:-06}"

dirname=fips-${loc}/${fdt}
mkdir -p $dirname
cd $dirname

dvc run \
    -w ../.. \
    -d fit-infection-kalman.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -d hopkins-vaccine/$fdt \
    -d healthdata/${fdt}/${loc} \
    -d covidcast-safegraph-home-prop-7dav/${fdt} \
    -o fits/${fdt}-fips${loc}/fit.RData \
    --force \
    -n fit-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript fit-infection-kalman.R
exit 0
dvc run \
    -w ../.. \
    -d write-formatted-forecast.R \
    -d covidhub-common.R \
    -d fits/${fdt}-fips${loc}/fit.RData \
    -o forecasts/${fdt}-fips${loc} \
    -o plots/${fdt}-fips${loc} \
    --force \
    -n forecast-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript write-formatted-forecast.R