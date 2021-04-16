#!/usr/bin/env bash

docker run \
	-d \
	--rm \
	-p 8788:8787 \
	-e USERID=1005 \
	-e PASSWORD=foo \
	--mount type=bind,src=$HOME/src/random-walks,dst=/home/rstudio/work \
	eamon/sir-kf:2021-04-16 /init
