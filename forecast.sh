#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-07-26}"
loc="${loc:-36}"

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