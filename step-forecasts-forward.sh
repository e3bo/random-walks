#!/usr/bin/bash

set -e

fdtstart=2021-02-08; 
for fdt in $(cat remaining-fdts); 
  do fdt=$fdt ./data-pull.sh && loc=06 fdt=$fdt fdtstart=$fdtstart ./forecast.sh; 
  fdtstart=$fdt; 
done