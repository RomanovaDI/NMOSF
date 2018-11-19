/*
matrix A and vector B structure in avalanche case:
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [0, 0, 0]
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [0, 1, 0]
.....................................................................
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [0, ny - 1, 0]
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [1, ny - 1, 0]
.....................................................................
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [nx - 1, ny - 1, 0]
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [nx - 1, ny - 1, 1]
.....................................................................
velocity_x, velocity_y, velocity_z, phase_fraction, pressure in cell [nx - 1, ny - 1, nz - 1]
*/

/*
matrix A and vector B structure in termogas case:
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [0, 0, 0]
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [0, 1, 0]
...............................................................................................................................................................................
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [0, ny - 1, 0]
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [1, ny - 1, 0]
...............................................................................................................................................................................
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [nx - 1, ny - 1, 0]
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [nx - 1, ny - 1, 1]
...............................................................................................................................................................................
concentration_N2, concentration_O2, concentration_CO2, concentration_H2O, pressure, saturation_water, saturation_oil, saturation_gas, temperature_flow, temperature_environment [nx - 1, ny - 1, nz - 1]
*/

/*
list of objects
   density_velocity
   density_velocity_velocity
   pressure
   shear_stress
   gravity_force
   grad_pressure
   div_density_velocity_velocity
   snow_volume_fraction
   snow_volume_fraction_velocity
*/

/*
list of numerical_schemes
	crank_nikolson
	half_forward_euler
	half_backward_euler
	forward_euler - unsteady solution
	backward_euler - steady solution
*/

/*
list of approximation orders
	first
	second
*/

/*
list of solution modes
	combined
	separated
	ultra_combined - don't working
*/

/*
list of solution methods
	VOF - volume of fluid method or finite volume method
	FDM - finite difference method
*/

/*
list of conditions types
	initial
	boundary
*/

/*
list of conditions objects
	velocity
	phase_fraction
	pressure
*/

/*
list of boundary condition modes
	zero_gradient_on_up_and_sides_no_slip_on_low
	zero_gradient_on_all
	fixed_value_on_up_and_sides_zero_gradient_on_low
	fixed value_on_all
*/

/*
list of initial conditions
	fixed_value
	ascii_map
	fixed_value_with_hydrostatic_pressure
	mass_is_mooving
*/

/*
numbering faces of cells is so
0 face - from beginning of x-axis
1 face - from end of x-axis
2 face - from beginning of y-axis
3 face - from end of y-axis
4 face - from beginning of z-axis
5 face - from end of z-axis
*/

#include "init_data.h"
#include "read_files.h"
#include "mesh_operations.h"
#include "utils.h"
#include "vtk_map_functions.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "array_functions.h"
#if AVALANCHE
#include "x_crank_nikolson_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_FDM_avalanche.h"
#include "x_backward_euler_second_combined_VOF_avalanche.h"
#include "x_backward_euler_second_combined_FDM_avalanche.h"
#include "t_second_combined_VOF_avalanche.h"
#endif
#if TERMOGAS
#include "t_second_combined_FDM_termogas.h"
#include "x_backward_euler_second_combined_FDM_termogas.h"
#include "t_second_separated_FDM_termogas.h"
#include "x_backward_euler_second_separated_FDM_termogas.h"
#endif
#include "t_test.h"
#include "create_matrix.h"
#include "matrix_functions.h"
#include "slu_ddefs.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

void display_usage(void);

/* modes of matrix creating */
#define PATTERN 0
#define MATRIX 1

#define SET_CONDITION(type, object, mode) SET_##type##_CONDITION_##object##_##mode(I)

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S, -t must be declared\n\n");
	return;
}

int main(int argc, char **argv)
{
	double t;
	time_t time1, time2;
	in II;
	in *I = &II;
	MPI_Init(&argc, &argv);
	if (set_parameters_termogas(I)) goto error;
	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime();
	if (read_asc_and_declare_variables(I)) goto error;
	if (I->my_rank == 0)
		if (do_interpolation(I)) goto error;
	if (I->nproc > 1)
		if (share_dx(I)) goto error;
	if (do_decomposition(I)) goto error;
#if AVALANCHE
	if (set_arrays(I)) goto error;
#endif
	if (make_boundary(I)) goto error;
	I->system_dimension = I->n_cells_multipl * I->nz * I->num_parameters;
	I->system_dimension_with_boundary = I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * I->num_parameters;
#if DEBUG
	printf("Process %d: I->system_dimension = %d, I->system_dimension_with_boundary = %d\n", I->my_rank, I->system_dimension, I->system_dimension_with_boundary);
#endif
	if ((I->B_prev = (double *) malloc(I->system_dimension_with_boundary * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->B_prev, 0, I->system_dimension_with_boundary * sizeof(double));
	if ((I->porousness = (double *) malloc(I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->porousness, 0, I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * sizeof(double));
	SET_CONDITION(initial, termogas, fixed_value);
	SET_CONDITION(boundary, termogas, no_boundaries_4_in_1_out);
	I->flag_first_time_step = 1;
	printf("Injection well: (%d, %d), (%d, %d), (%d, %d), (%d, %d).\n", I->nx / 2, I->ny / 2, I->nx / 2, I->ny / 2 - 1, I->nx / 2 - 1, I->ny / 2, I->nx / 2 - 1, I->ny / 2 - 1);
	I->marker_of_step = 0;
	while (I->time < I->end_time) {
		if (I->my_rank == 0) {
			printf("Time is %lf of %lf, time step %d\n", I->time, I->end_time, I->time_step);
		}
		if ((I->nproc > 1) && (reconstruct_src(I))) return 1;
		if (set_array_of_parameters_termogas(I)) return 1;
		if (check_convection_speed(I)) return 1;
		if (I->nan_flag) {
			printf("There are NANs\n");
			goto error;
		}
		if (I->negative_num_flag) {
			printf("There are negative numbers in values from set of positive numbers\n");
			goto error;
		}
		if ((I->flag_first_time_step) && (I->nproc > 1) && (reconstruct_src(I))) return 1;
#if DEBUG
		if ((I->my_rank == 0) && (I->nproc > 1) && (print_gl_B(I, 4))) return 1;
#endif
		time(&time1);
		if ((I->marker_of_step % 2) || (I->flag_first_time_step)) {
			if (print_vtk(I)) {
				printf("Error printing vtk file\n");
				goto error;
			}
		}
		time(&time2);
		printf("Time of printing data is %lfsec.\n", difftime(time2, time1));
//		if (i % 100 == 0) {
//			if (print_vtk(I, i / 100) == 1) {
//				printf("Error printing vtk file\n");
//				goto error;
//			}
//		}
#if DEBUG
		if (print_parameter_in_subdomains(I, 4)) return 1;
#endif
#if AVALANCHE
		if (create_Ab_avalanche(I) == 1) goto error;
#endif
#if TERMOGAS
		//if (!I->flag_first_time_step)
		//	calculate_disparity_pressure(I);
		I->courant_number = 0.7;
		I->avarage_velocity_global_value = avarage_velocity_global(I) / I->courant_number;
		I->dt = I->dx[0] / I->avarage_velocity_global_value;
		printf("avarage_velocity_global = %lf\tdt = %lf\n", I->avarage_velocity_global_value, I->dt);
		if (I->flag_first_time_step)
			if (create_Ab_termogas(I)) goto error;
		time(&time1);
		int tmp = 0;
		do {
			printf("Matrix creating iteration %d\n", tmp);
			I->courant_flag = 0;
			if (fill_Ab_termogas(I)) goto error;
			if (I->courant_flag) {
				I->avarage_velocity_global_value *= 2;
				I->dt = I->dx[0] / I->avarage_velocity_global_value;
			}
			tmp++;
		} while (I->courant_flag);
		time(&time2);
		printf("Time of matrix creating is %lfsec.\n", difftime(time2, time1));
#endif
		//print_A_csr(I);
		//printf("1\n");
		//check_values_in_B(I);
		time(&time1);
		if (solve_matrix(I)) goto error;
		time(&time2);
		//printf("2\n");
		check_values_in_B(I);
		printf("Time of matrix solving is %lfsec.\n", difftime(time2, time1));
		if (I->my_rank == 0) {
			if (print_oil_production(I)) goto error;
		}
		//if (write_B_to_B_prev(I)) goto error;
		if (I->marker_of_step % 2 == 0) {
			if (write_pressure_and_temperature(I)) goto error;
			I->marker_of_step++;
		} else {
			if (write_B_to_B_prev(I)) goto error;
			I->time += I->dt;
			I->time_step++;
			I->marker_of_step++;
		}
		if (I->flag_first_time_step)
			I->flag_first_time_step = 0;
		//printf("3\n");
		//check_values_in_B_prev(I);
		//printf("1\n");
		//check_for_symmetry(I);
		if (check_sum(I)) goto error;
		if (I->nan_flag) {
			printf("There are NANs\n");
			goto error;
		}
		SET_CONDITION(boundary, termogas, no_boundaries_4_in_1_out);
		//printf("2\n");
		//check_for_symmetry(I);
		//printf("4\n");
		//check_values_in_B_prev(I);
	}
	time(&time1);
	if (print_vtk(I)) {
		printf("Error printing vtk file\n");
		goto error;
	}
	time(&time2);
	save_I_termogas(I);
	if (free_massives(I)) goto error;

	if (I->my_rank == 0) printf("Calculations finished successfully\n");
	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
	printf("Processor %d PID %d: work_time = %lf sec\n", I->my_rank, getpid(), t);
	MPI_Finalize();
	return 0;
error:
	printf("Error in process %d\n", I->my_rank);
//	display_usage();
	free_massives(I);
	MPI_Finalize();
	return 1;
}
//if (solve_test_matrix()) goto error;
//if (set_parameters_avalanche(I)) goto error;
//SET_CONDITION(initial, phase_fraction, fixed_value);
//SET_CONDITION(initial, pressure, fixed_value);
//SET_CONDITION(initial, velocity, fixed_value);
//SET_CONDITION(boundary, velocity, zero_gradient_on_y_and_x_and_upper_sides_no_slip_on_low); //mass
//SET_CONDITION(boundary, velocity, zero_gradient_on_y_sides_no_slip_on_other_upper_wall_is_mooving); //cavity
//SET_CONDITION(boundary, velocity, zero_gradient_on_y_and_x_sides_no_slip_on_other_upper_wall_is_mooving); //couette
//SET_CONDITION(boundary, phase_fraction, zero_gradient_on_all);
//SET_CONDITION(boundary, pressure, zero_gradient_on_all);
//if (barotropy_density(I)) goto error;
//if (barotropy_pressure(I)) goto error;
//if (check_conservation_of_mass(I)) {
//	printf("Mass conservation equation failed\n");
//}
