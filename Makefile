NAME	=	NMOSF
CC		=	mpicc
LD		=	mpicc
DEBUG	=	-g -O0
CFLAGS	=	$(DEBUG) -shared -fPIC -I src/lib/ -L build/lib/ -I /usr/include/superlu/ -L /usr/lib/x86_64-linux-gnu/
LDFLAGS	=	$(DEBUG) -I src/lib/ -L build/lib/ -I /usr/include/superlu/ -L /usr/lib/x86_64-linux-gnu/

OBJ0	=	init_data
LIB0	=	-lm

OBJ1	=	read_files utils
LIB1	=	$(LIB0) $(addprefix -l,$(OBJ0))

OBJ2	=	mesh_operations vtk_map_functions boundary_conditions initial_conditions array_functions t_test t_second_combined_FDM_termogas x_backward_euler_second_combined_FDM_termogas t_second_separated_FDM_termogas x_backward_euler_second_separated_FDM_termogas
LIB2	=	$(LIB1) $(addprefix -l,$(OBJ1))

OBJ3	=	matrix_functions
LIB3	=	$(LIB2) $(addprefix -l,$(OBJ2)) -lsuperlu

OBJ4	=	create_matrix
LIB4	=	$(LIB3) $(addprefix -l,$(OBJ3))

LDLIB	=	$(LIB4) $(addprefix -l,$(OBJ4))

#export NMOSF=/home/daria/cases/NMOSF/

all: $(NAME)

$(NAME): 
	$(CC) $(CFLAGS) src/lib/init_data.c $(LIB0) -o build/lib/libinit_data.so
	$(CC) $(CFLAGS) src/lib/read_files.c $(LIB1) -o build/lib/libread_files.so
	$(CC) $(CFLAGS) src/lib/utils.c $(LIB1) -o build/lib/libutils.so
	$(CC) $(CFLAGS) src/lib/mesh_operations.c $(LIB2) -o build/lib/libmesh_operations.so
	$(CC) $(CFLAGS) src/lib/vtk_map_functions.c $(LIB2) -o build/lib/libvtk_map_functions.so
	$(CC) $(CFLAGS) src/lib/boundary_conditions.c $(LIB2) -o build/lib/libboundary_conditions.so
	$(CC) $(CFLAGS) src/lib/initial_conditions.c $(LIB2) -o build/lib/libinitial_conditions.so
	$(CC) $(CFLAGS) src/lib/array_functions.c $(LIB2) -o build/lib/libarray_functions.so
	#gcc -g -O0 -shared src/lib/t_second_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libt_second_combined_VOF_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/t_second_ultra_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libt_second_ultra_combined_VOF_avalanche.so -fPIC
	$(CC) $(CFLAGS) src/lib/t_test.c $(LIB2) -o build/lib/libt_test.so
	#gcc -g -O0 -shared src/lib/x_crank_nikolson_second_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_crank_nikolson_second_combined_VOF_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/x_forward_euler_second_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_forward_euler_second_combined_VOF_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/x_forward_euler_second_combined_FDM_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_forward_euler_second_combined_FDM_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/x_backward_euler_second_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_backward_euler_second_combined_VOF_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/x_backward_euler_second_combined_FDM_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_backward_euler_second_combined_FDM_avalanche.so -fPIC
	#gcc -g -O0 -shared src/lib/x_backward_euler_second_ultra_combined_VOF_avalanche.c -lm -I src/lib/ -L build/lib/ -linit_data -lutils -o build/lib/libx_backward_euler_second_ultra_combined_VOF_avalanche.so -fPIC
	$(CC) $(CFLAGS) src/lib/t_second_combined_FDM_termogas.c $(LIB2) -o build/lib/libt_second_combined_FDM_termogas.so
	$(CC) $(CFLAGS) src/lib/x_backward_euler_second_combined_FDM_termogas.c $(LIB2) -o build/lib/libx_backward_euler_second_combined_FDM_termogas.so
	$(CC) $(CFLAGS) src/lib/t_second_separated_FDM_termogas.c $(LIB2) -o build/lib/libt_second_separated_FDM_termogas.so
	$(CC) $(CFLAGS) src/lib/x_backward_euler_second_separated_FDM_termogas.c $(LIB2) -o build/lib/libx_backward_euler_second_separated_FDM_termogas.so
	$(CC) $(CFLAGS) src/lib/matrix_functions.c $(LIB3) -o build/lib/libmatrix_functions.so
	$(CC) $(CFLAGS) src/lib/create_matrix.c $(LIB4) -o build/lib/libcreate_matrix.so
	$(LD) $(LDFLAGS) src/asc2Ab.c $(LDLIB) -o build/NMOSF
	export LD_LIBRARY_PATH=./build/lib
#./build/NMOSF
#valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes ./build/NMOSF
#gdb ./build/NMOSF
#rm tmp/*
#gcc -g -O0 asc2Ab.c -lm -I/usr/include/superlu/ -L/usr/lib/x86_64-linux-gnu/ -lsuperlu
#valgrind --leak-check=full --leak-resolution=med --show-leak-kinds=all --track-origins=yes
#./a.out -m map_for_verification_Ab_1.asc -r map_for_verification_Ab_region_1.asc -H 15 -D 7 -x 1 -y 1 -z 1 -a 200 -s 1 -p 101325 -v 0.0000148 -k 1000 -i 0.5 -l 2000 -S 0.00001 -t 0.3


clean:
	rm -rf result
	mkdir result
	rm -rf build
	mkdir build
	mkdir build/lib
#	mkdir tmp

#	-lt_second_combined_VOF_avalanche \
#	-lt_second_ultra_combined_VOF_avalanche \
#	-lx_crank_nikolson_second_combined_VOF_avalanche \
#	-lx_forward_euler_second_combined_VOF_avalanche \
#	-lx_forward_euler_second_combined_FDM_avalanche \
#	-lx_backward_euler_second_combined_VOF_avalanche \
#	-lx_backward_euler_second_combined_FDM_avalanche \
#	-lx_backward_euler_second_ultra_combined_VOF_avalanche \
