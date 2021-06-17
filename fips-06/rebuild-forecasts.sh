#!/usr/bin/bash

set -e

for fdt in $(cat wave-fdts); do 
  cd $fdt
  dvc repro -s forecast-${fdt}-06; 
  cd ..
done
