#!/usr/bin/env bash

docker run \
	-d \
	--rm \
	-p 8788:8787 \
	-e USERID=1001 \
	-e PASSWORD=foo \
	--mount type=bind,src=$HOME/src/random-walks,dst=/home/rstudio/work \
        --mount type=bind,src=$HOME/src/covid19-forecast-hub,dst=/home/rstudio/hub \
	docker.io/eamon/sir-kf:2021-06-24 /init
