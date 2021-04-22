#!/usr/bin/env bash

set -e


fdtstart="${fdtstart:-2020-06-29}"
fdt="${fdt:-2020-06-29}"
loc="${loc:-36}"

dirname=fips-${loc}/${fdt}
mkdir -p $dirname
cd $dirname

if [$fdtstart == $fdt]; then
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
      fdtstart=$fdtstart fdt=$fdt loc=$loc Rscript fit-infection-kalman.R
else 
  dvc run \
      -w ../.. \
      -d fit-infection-kalman.R \
      -d covidhub-common.R \
      -d hopkins/$fdt \
      -d covidcast-safegraph-home-prop/${fdt} \
      -d fits/${fdtstart}-fips${loc}/fit.RData \
      -o fits/${fdt}-fips${loc}/fit.RData \
      -M fits/${fdt}-fips${loc}/fit-metrics.json \
      --force \
      -n fit-$fdt-$loc \
      fdtstart=$fdtstart fdt=$fdt loc=$loc Rscript fit-infection-kalman.R
fi

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