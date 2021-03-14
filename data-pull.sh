#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"
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