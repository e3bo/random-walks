#!/usr/bin/env bash

set -e

fdt="${fdt:-2021-03-29}"
locs=("06")

dvc run \
    -d pull-hopkins-ts-from-date.R \
    -o hopkins/$fdt \
    -n hopkins-$fdt \
    --force \
    fdt=$fdt ./pull-hopkins-ts-from-date.R

if [[ ${fdt} < "2020-11-16" ]]; then
    exit 0
fi

#locs=("01" "02" "04" "05" "06" "08" "09" "10" "11" "12" "13" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "44" "45" "46" "47" "48" "49" "50" "51" "53" "54" "55" "56")
for loc in "${locs[@]}"
do
dvc run \
    -d pull-healthdata-ts-from-date.R \
    -o healthdata/$fdt/$loc \
    -n healthdata-$fdt-$loc \
    --force \
    fdt=$fdt loc=$loc ./pull-healthdata-ts-from-date.R
done

if [[ ${fdt} < "2021-02-12" ]]; then
    exit 0
fi

dvc run \
    -d pull-hopkins-vaccine-ts-from-date.R \
    -o hopkins-vaccine/$fdt \
    -n hopkins-vaccine-$fdt \
    --force \
    fdt=$fdt ./pull-hopkins-vaccine-ts-from-date.R