[![DOI](https://zenodo.org/badge/285368880.svg)](https://zenodo.org/badge/latestdoi/285368880)

This repository contains the code of a forecaster for [the COVID 19 forecast hub](https://covid19forecasthub.org/), which collects forecasts of COVID 19 indicators,
such as official reports of deaths and cases, in a standardized format. In brief, this forecaster is a state space SEIRD model which allows for parameters to change over time according to a random walk.

A workflow for running and evaluating this forecaster with versioned data has been encoded with the [DVC](https://dvc.org) system. It may be run on a system built according to the Dockerfile in the `docker` subdirectory using the following command line: `dvc repro analyze-scores`. A manuscript introducing this forecaster and presenting the evaluation should be available soon.

The above command `dvc repro analyze-scores` works by reading instructions on how to create the resulsts from the file dvc.yaml. 
In this file, the entire workflow is specified in a directed hierarchy of key-value pairs.
At the top level is the 'stage' key, and at the next level is the name of each stage in the workflow.
Each stage comprises the keys 'cmd', 'deps', and 'outs'.
The value of 'cmd' for a stage is the command line used to run the stage.
The value of 'deps' for a stage is a list of dependences for a stage.
Any changes in these files or directories necessitate a re-running of the command line to ensure that the output of the stage is up-to-date.
The value of 'outs' is a list of the files and directories produced by the stage.
Thefore, even without the dvc program, it is possible by examining the file `dvc.yaml` to determine the sequence of commands needed to run the 'analyze-scores' stage which produces the main results of the paper.

However, the number of stages and files involved is large and so we next offer an overview of the organization of the code that may help users find a particular item of interest.
The  log likelihood functions are contained in the file `InfectinKalman.jl`.
These functions are called from the script `fit-infection-kalman.R` to fit the model to a particular version of data and calculate simple simple metrics.
Note that an R version of the log likelihood function is used to calculate some of the metrics.
The R and Julia functions should be identical for the outputs which they both produce, while the R version produces additional outputs that are not needed during the optimization.
The creation of forecasts in the COVID-19 Forecast Hub format is done by `write-formatted-forecast.R`, which also contains code for making a table of parameter estimates and many plots of the forecasts.
Calculation of weighted interval scores is done in `analyze-scores.R`.
See the `data-pull.sh` bash script for an example for how to pull the cases, deaths, and hosptalization data for a particular week.
The URLs of the Google Community Mobility reports used are obtained with `lookup-wayback-urls.R`, then the files are downloaded and filtered to US state data (to save space) by `download-google-mob.R`.
