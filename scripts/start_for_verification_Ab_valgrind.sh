#!/bin/bash
gcc -g -O0 asc2Ab.c -lm
valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes ./a.out -m map_for_verification_Ab_1.asc -r map_for_verification_Ab_region_1.asc -H 15 -D 5 -x 1 -y 1 -z 1 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001
