#!/bin/sh
set -x

# run the built container, mount current folder into the container
docker run --rm -it -v $PWD:/home/app/voxcraft-sim voxcraft-sim /bin/sh
