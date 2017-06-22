/*
matrix A and vector B structure:
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
	forward_euler
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
#include "x_crank_nikolson_second_combined_VOF.h"
#include "x_forward_euler_second_combined_VOF.h"
#include "x_forward_euler_second_combined_FDM.h"
#include "t_second_combined_VOF.h"
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
	int i, time_steps;
	in II;
	in *I = &II;
	//if (solve_test_matrix()) goto error;
	if (set_parameters(I)) goto error;
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
	SET_CONDITION(initial, phase_fraction, fixed_value);
	SET_CONDITION(initial, pressure, fixed_value);
	SET_CONDITION(initial, velocity, fixed_value);
	time_steps = I->end_time / I->dt;
	I->flag_first_time_step = 1;
	for (i = 0; i <= time_steps; i++) {
		SET_CONDITION(boundary, velocity, zero_gradient_on_y_sides_no_slip_on_other_upper_wall_is_mooving);
		SET_CONDITION(boundary, phase_fraction, zero_gradient_on_all);
		SET_CONDITION(boundary, pressure, zero_gradient_on_all);
		if (print_vtk(I, i) == 1) {
			printf("Error printing vtk file\n");
			goto error;
		} else {
			printf("Result printed to vtk file\n");
		}
		if (create_Ab(I) == 1) goto error;
		if (i == 0)
			I->flag_first_time_step = 0;
		if (solve_matrix(I)) goto error;
		if (barotropy_pressure(I)) goto error;
		write_B_to_B_prev(I);
		if (check_conservation_of_mass(I)) {
			printf("Mass conservation equation failed\n");
		}
		if (i == time_steps) {
			if (print_vtk(I, i + 1) == 1) {
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
