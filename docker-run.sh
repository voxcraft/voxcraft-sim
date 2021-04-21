#!/bin/sh
set -x

# run the built container, mount current folder into the container
docker run --rm -it --gpus all -v $PWD:/home/app/voxcraft-sim voxcraft-sim /bin/bash
