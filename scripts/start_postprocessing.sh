#!/bin/bash
set -x
export LD_LIBRARY_PATH=./build/lib
mpirun -np 1 ./build/postprocessing 2>&1 | tee result/postprocessing.log
