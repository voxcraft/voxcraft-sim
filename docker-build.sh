#!/bin/sh
set -x

# build a container based on the Dockerfile
# pass in the current uid and gid, so the user inside the docker will align with the current system
docker build --build-arg UID=`id -u` --build-arg GID=`id -g` -t voxcraft-sim .
