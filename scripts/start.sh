#!/bin/bash
set -x
rm -rf result
mkdir result
mkdir tmp
mkdir build
mkdir build/lib
gcc -g -O0 -shared src/lib/init_data.c -lm -o build/lib/libinit_data.so -fPIC
gcc -g -O0 -shared src/lib/read_files.c -lm -o build/lib/libread_files.so -fPIC
gcc -g -O0 -shared src/lib/mesh_operations.c -lm -o build/lib/libmesh_operations.so -fPIC
gcc -g -O0 -shared src/lib/utils.c -lm -I/usr/include/superlu/ -L/usr/lib/x86_64-linux-gnu/ -lsuperlu -o build/lib/libutils.so -fPIC
gcc -g -O0 -shared src/lib/vtk_map_functions.c -lm -o build/lib/libvtk_map_functions.so -fPIC
gcc -g -O0 -shared src/lib/boundary_conditions.c -lm -o build/lib/libboundary_conditions.so -fPIC
gcc -g -O0 -shared src/lib/initial_conditions.c -lm -o build/lib/libinitial_conditions.so -fPIC
gcc -g -O0 -shared src/lib/array_functions.c -lm -o build/lib/libarray_functions.so -fPIC
gcc -g -O0 -shared src/lib/t_second_combined_VOF_avalanche.c -lm -o build/lib/libt_second_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/t_second_ultra_combined_VOF_avalanche.c -lm -o build/lib/libt_second_ultra_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/t_test.c -lm -o build/lib/libt_test.so -fPIC
gcc -g -O0 -shared src/lib/x_crank_nikolson_second_combined_VOF_avalanche.c -lm -o build/lib/libx_crank_nikolson_second_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/x_forward_euler_second_combined_VOF_avalanche.c -lm -o build/lib/libx_forward_euler_second_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/x_forward_euler_second_combined_FDM_avalanche.c -lm -o build/lib/libx_forward_euler_second_combined_FDM_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/x_backward_euler_second_combined_VOF_avalanche.c -lm -o build/lib/libx_backward_euler_second_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/x_backward_euler_second_combined_FDM_avalanche.c -lm -o build/lib/libx_backward_euler_second_combined_FDM_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/x_backward_euler_second_ultra_combined_VOF_avalanche.c -lm -o build/lib/libx_backward_euler_second_ultra_combined_VOF_avalanche.so -fPIC
gcc -g -O0 -shared src/lib/create_matrix.c \
	-I src/lib/ \
	-L build/lib \
	-lm \
	-lt_second_combined_VOF_avalanche \
	-lt_second_ultra_combined_VOF_avalanche \
	-lx_crank_nikolson_second_combined_VOF_avalanche \
	-lx_forward_euler_second_combined_VOF_avalanche \
	-lx_forward_euler_second_combined_FDM_avalanche \
	-lx_backward_euler_second_combined_VOF_avalanche \
	-lx_backward_euler_second_combined_FDM_avalanche \
	-lx_backward_euler_second_ultra_combined_VOF_avalanche \
	-lt_test \
	-o build/lib/libcreate_matrix.so \
	-fPIC
gcc -g -O0 -shared src/lib/matrix_functions.c -lm -I /usr/include/superlu/ -L /usr/lib/x86_64-linux-gnu/ -lsuperlu -o build/lib/libmatrix_functions.so -fPIC
gcc -g -O0  src/asc2Ab.c \
	-I /usr/include/superlu/ -I src/lib/ \
	-L /usr/lib/x86_64-linux-gnu/ -L build/lib/ \
	-lm \
	-linit_data \
	-lread_files \
	-lmesh_operations \
	-lutils \
	-lvtk_map_functions \
	-lboundary_conditions \
	-linitial_conditions \
	-larray_functions \
	-lt_second_combined_VOF_avalanche \
	-lt_second_ultra_combined_VOF_avalanche \
	-lt_test \
	-lx_crank_nikolson_second_combined_VOF_avalanche \
	-lx_forward_euler_second_combined_VOF_avalanche \
	-lx_forward_euler_second_combined_FDM_avalanche \
	-lx_backward_euler_second_combined_VOF_avalanche \
	-lx_backward_euler_second_combined_FDM_avalanche \
	-lx_backward_euler_second_ultra_combined_VOF_avalanche \
	-lcreate_matrix \
	-lmatrix_functions \
	-o build/NMOSF
export LD_LIBRARY_PATH=./build/lib
./build/NMOSF
#valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes ./build/NMOSF
#gdb ./build/NMOSF
#rm tmp/*
#gcc -g -O0 asc2Ab.c -lm -I/usr/include/superlu/ -L/usr/lib/x86_64-linux-gnu/ -lsuperlu
#valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes
#./a.out -m map_for_verification_Ab_1.asc -r map_for_verification_Ab_region_1.asc -H 15 -D 7 -x 1 -y 1 -z 1 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001 -t 0.3
