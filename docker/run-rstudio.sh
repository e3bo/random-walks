#!/usr/bin/env bash

podman run \
	-d \
	--rm \
	-p 8788:8787 \
	-e USERID=$(id -u) \
	-e PASSWORD=foo \
	--mount type=bind,src=$HOME/src/random-walks,dst=/home/rstudio/work \
        --mount type=bind,src=$HOME/src/covid19-forecast-hub,dst=/home/rstudio/hub \
	docker.io/eamon/sir-kf:2021-07-17 /init
