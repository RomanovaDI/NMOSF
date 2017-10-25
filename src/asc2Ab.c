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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

struct timeval tv1,tv2,dtv;
struct timezone tz;

void display_usage(void);
void time_start(void);
long time_stop(void);

/* modes of matrix creating */
#define PATTERN 0
#define MATRIX 1

#define SET_CONDITION(type, object, mode) SET_##type##_CONDITION_##object##_##mode(I)

void time_start(void)
{
	gettimeofday(&tv1, &tz);
	return;
}

long time_stop(void)
{
	gettimeofday(&tv2, &tz);
	dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
	dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
	if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
		return dtv.tv_sec*1000+dtv.tv_usec/1000;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S, -t must be declared\n\n");
}

int main(int argc, char **argv)
{
	int i, time_steps, j;
	in II;
	in *I = &II;
//	solve_test_matrix();
//	return 0;
	if (set_parameters_termogas(I)) goto error;
	if (read_asc_and_declare_variables(I)) goto error;
	if (do_interpolation(I)) goto error;
	if (set_arrays(I)) goto error;
	if (make_boundary(I)) goto error;
	I->system_dimension = I->n_cells_multipl * I->nz * I->num_parameters;
	I->system_dimension_with_boundary = I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * I->num_parameters;
	if ((I->B_prev = (double *) malloc(I->system_dimension_with_boundary * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->B_prev, 0, I->system_dimension_with_boundary * sizeof(double));
	SET_CONDITION(initial, termogas, fixed_value);
	time_steps = I->end_time / I->dt;
	I->flag_first_time_step = 1;
//	time_steps = 20;
	for (i = 0; i <= time_steps; i++) {
		printf("Time step %d of %d\n", i + 1, time_steps);
		I->time_step = i;
		SET_CONDITION(boundary, termogas, no_bounadries_4_in_1_out);
		if (i % 10 == 0) {
			if (print_vtk(I, i / 10) == 1) {
				printf("Error printing vtk file\n");
				goto error;
			} else {
				printf("Result printed to vtk file\n");
			}
		}
		//for (j = 0; j < 5; j++) {
		//for (j = 0; j < 2; j++) {
		//	printf("Equation %d\n", j);
		//	I->equation_num = j;
#if AVALANCHE
			if (create_Ab_avalanche(I) == 1) goto error;
#endif
#if TERMOGAS
			if (create_Ab_termogas(I) == 1) goto error;
#endif
			if (i == 0)
				I->flag_first_time_step = 0;
			if (solve_matrix(I)) goto error;
			if (print_oil_production(I)) goto error;
			if (write_B_to_B_prev(I)) goto error;
			if (check_sum(I)) goto error;
//			if (print_vtk(I, j + (i + 1) * 1000) == 1) {
//				printf("Error printing vtk file\n");
//				goto error;
//			} else {
//				printf("Result printed to vtk file\n");
//			}
		//}
		if (i == time_steps) {
			if (print_vtk(I, i / 10 + 1) == 1) {
				printf("Error printing vtk file\n");
				goto error;
			} else {
				printf("Result printed to vtk file\n");
			}
		}
	}
	if (free_massives(I) == 1) goto error;

	printf("Calculations finished successfully\n");
	return 0;
error:
	printf("Error\n");
	display_usage();
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
