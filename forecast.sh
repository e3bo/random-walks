#!/usr/bin/env bash

set -e

fdt="${fdt:-2021-03-29}"
loc="${loc:-36}"

dirname=fips-${loc}/${fdt}
mkdir -p $dirname
cd $dirname

dvc run \
    -w ../.. \
    -d fit-infection-kalman.R \
    -d covidhub-common.R \
    -d hopkins/$fdt \
    -d covidcast-safegraph-home-prop/${fdt} \
    -o fits/${fdt}-fips${loc}/fit.RData \
    -M fits/${fdt}-fips${loc}/fit-metrics.json \
    --force \
    -n fit-$fdt-$loc \
    fdt=$fdt loc=$loc Rscript fit-infection-kalman.R

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