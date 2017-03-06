#!/bin/bash
gcc -g asc2Ab.c -lm
gdb ./a.out -m map_for_verification_Ab.asc -r map_for_verification_Ab_region.asc -H 10 -D 2 -x 2 -y 2 -z 2 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001
