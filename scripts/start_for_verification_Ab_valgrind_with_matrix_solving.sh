#!/bin/bash
set -x
gcc -shared src/lib/init_data.c -lm -o build/lib/libinit_data.so -fPIC
gcc -shared src/lib/read_files.c -lm -o build/lib/libread_files.so -fPIC
gcc -shared src/lib/mesh_operations.c -lm -o build/lib/libmesh_operations.so -fPIC
gcc -shared src/lib/utils.c -lm -I/usr/include/superlu/:/home/daria/cases/NMOSF/src/lib -L/usr/lib/x86_64-linux-gnu/:/home/daria/cases/NMOSF/src/lib -lsuperlu -o build/lib/libutils.so -fPIC
gcc -shared src/lib/boundary_conditions.c -lm -o build/lib/libboundary_conditions.so -fPIC
gcc -shared src/lib/initial_conditions.c -lm -o build/lib/libinitial_conditions.so -fPIC
gcc -shared src/lib/array_functions.c -lm -o build/lib/libarray_functions.so -fPIC
gcc -shared src/lib/t_second_combined_VOF.c -lm -o build/lib/libt_second_combined_VOF.so -fPIC
gcc -shared src/lib/x_crank_nikolson_second_combined_VOF.c -lm -o build/lib/libx_crank_nikolson_second_combined_VOF.so -fPIC
#gcc -g -O0 asc2Ab.c -lm -I/usr/include/superlu/ -L/usr/lib/x86_64-linux-gnu/ -lsuperlu
#valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes
#./a.out -m map_for_verification_Ab_1.asc -r map_for_verification_Ab_region_1.asc -H 15 -D 7 -x 1 -y 1 -z 1 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001 -t 0.3
