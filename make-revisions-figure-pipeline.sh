#!/bin/env bash

dvc run -d hopkins -d study-indicators-revisions.R -d covidhub-common.R -o indicators-revisions-cal.png --name make-revisions-plot ./study-indicators-revisions.R
