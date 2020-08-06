#!/usr/bin/env bash

docker run \
	-d \
	--rm \
	-p 8788:8787 \
	-e USERID=1005 \
	-e DISABLE_AUTH=true \
	--mount type=bind,src=$HOME/src/random-walks,dst=/home/rstudio/work \
	eamon/random-walks:2020-08-06 /init
