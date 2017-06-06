#!/bin/bash
gcc asc2Ab.c -lm
./a.out -m test_map.asc -r test_map_snow_cover.asc -H 20 -D 2 -x 2 -y 2 -z 2 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001
