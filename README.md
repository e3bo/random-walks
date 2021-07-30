[![DOI](https://zenodo.org/badge/285368880.svg)](https://zenodo.org/badge/latestdoi/285368880)

This repository contains the code of a forecaster for [the COVID 19 forecast hub](https://covid19forecasthub.org/), which collects forecasts of COVID 19 indicators,
such as official reports of deaths and cases, in a standardized format. In brief, this forecaster is a state space SEIRD model which allows for parameters to change over time according to a random walk.

A workflow for running and evaluating this forecaster with versioned data has been encoded with the [DVC](https://dvc.org) system. It may be run on a system built according to the Dockerfile in the `docker` subdirectory using the following command line: `dvc repro analyze-scores`. A manuscript introducing this forecaster and presenting the evaluation should be available soon.
