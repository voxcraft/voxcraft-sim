#!/bin/sh
cwd=$(pwd)
BuildSim=true
RebuildAll=false
while getopts “:f” opt; do
  case $opt in #-f means rebuild from scratch
    f) RebuildAll=true ;;
  esac
done


if $BuildSim; then
    # For voxcraft-sim (Server side)
    cd ../..
    if $RebuildAll; then
        echo "Rebuilding voxcraft-sim."
        rm build_Release/ -rf
        mkdir build_Release
    else
        echo "Making voxcraft-sim."
    fi
    cd build_Release
    if $RebuildAll; then
        cmake -DCMAKE_BUILD_TYPE=Release -DCUDA_DEBUG=OFF ..
    fi
    cmake --build . -j 7
    status=$?
    cd $cwd
    cp ../../build_Release/voxcraft-sim .
    cp ../../build_Release/vx3_node_worker .
    exit $status
fi