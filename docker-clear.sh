#!/bin/bash

docker kill $(docker ps -q)
docker rm $(docker ps -aq)
sudo service docker restart