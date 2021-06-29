#!/bin/env bash

dvc run -d fits \
-d compare-params-by-forecast-date.R \
-d covidhub-common.R \
-o res-effect-by-forecast-date.png \
--name beta-res-by-fdt \
./compare-params-by-forecast-date.R