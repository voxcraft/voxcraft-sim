#!/bin/sh
set -x

rm build -rf && mkdir build && cd build && cmake .. && make -j 7