#!/usr/bin/env bash

docker run \
	-d \
	--rm \
	-p 8788:8787 \
	-e USERID=1000 \
	-e DISABLE_AUTH=true \
	--mount type=bind,src=$HOME/src/random-walks,dst=/home/rstudio/work \
	eamon/sir-ekf:2021-01-08 /init
