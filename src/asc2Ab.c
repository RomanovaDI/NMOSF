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
   density_snow_volume_fraction
   density_snow_volume_fraction_velocity
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
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "x_crank_nikolson_second_comined_VOF.h"
#include "t_second_combine_VOF.h"
#include "slu_ddefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

struct timeval tv1,tv2,dtv;
struct timezone tz;

int print_vtk(int n);
void display_usage(void);
void time_start(void);
long time_stop(void);

/* modes of matrix creating */
#define PATTERN 0
#define MATRIX 1

#define SET_CONDITION(type, object, mode) SET_##type##_CONDITION_##object##_##mode()

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

int print_vtk(int n)
{
	int nn = n;
	int i = 0, j, k, a;
	while (nn > 0) {
		nn /= 10;
		i++;
	}
	//char file_name[8 + i];
	char file_name[20 + i];
	sprintf(file_name, "map%d.vtk", n);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "slope\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(f, "POINTS %d double\n", n_points_multipl * ((int) (hight / cellsize) * kz + 1));
	for (k = 0; k <= (int) (hight / cellsize) * kz; k++) {
		for (i = 0; i < (nrows - 1) * kx + 1; i ++) {
			for (j = 0; j < (ncols - 1) * ky + 1; j++) {
				if (ind_multipl[i * ((ncols - 1) * ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * cellsize / (double) kx, j * cellsize / (double) ky,
						mass1[i * ((ncols - 1) * ky + 1) + j] + k * cellsize / (double) kz);
			}
		}
	}
	
	a = 0;
	fprintf(f, "CELLS %d %d\n", n_cells * kx * ky * kz * (int) (hight / cellsize),
		n_cells * kx * ky * kz * (int) (hight / cellsize) * 9);
	for (k = 0; k < (int) (hight / cellsize) * kz; k++) {
		for (i = 1; i < (nrows - 1) * kx + 1; i++) {
			for (j = 1; j < (ncols - 1) * ky + 1; j++) {
				if ((ind_multipl[i * ((ncols - 1) * ky + 1) + j] != -1) &&
					(ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j] != -1) &&
					(ind_multipl[i * ((ncols - 1) * ky + 1) + j - 1] != -1) &&
					(ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j - 1] != -1)) {
						fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8,
							k * n_points_multipl + ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j - 1],
							k * n_points_multipl + ind_multipl[i * ((ncols - 1) * ky + 1) + j - 1],
							k * n_points_multipl + ind_multipl[i * ((ncols - 1) * ky + 1) + j],
							k * n_points_multipl + ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[i * ((ncols - 1) * ky + 1) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[i * ((ncols - 1) * ky + 1) + j],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j]);
				}
			}
		}
	}
	fprintf(f, "CELL_TYPES %d\n", n_cells * kx * ky * kz * (int) (hight / cellsize));
	for (i = 0; i < n_cells * kx * ky * kz * (int) (hight / cellsize); i++) {
		fprintf(f, "%d\n", 12);
	}

	fprintf(f, "CELL_DATA %d\n", n_cells_multipl * nz);
	fprintf(f, "VECTORS velocity double\n");
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1)
					fprintf(f, "%20.10f\t%20.10f\t%20.10f\n", B_prev[B_IND(0, i, j, k)], B_prev[B_IND(1, i, j, k)], B_prev[B_IND(2, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1)
					fprintf(f, "%20.10f\n", B_prev[B_IND(4, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS snow_volume_fraction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1)
					fprintf(f, "%20.10f\n", B_prev[B_IND(3, i, j, k)]);
			}
		}
	}
	fclose(f);
	sprintf(file_name, "B_prev%d.txt", n);
	f = fopen(file_name, "w");
	int p;
	for (p = 0; p < num_parameters; p++) {
		fprintf(f, "parameter %d\n", p);
		for (k = -stencil_size; k < nz + stencil_size; k++) {
			fprintf(f, "k = %d\n", k);
			for (i = -stencil_size; i < nx + stencil_size; i++) {
				for (j = -stencil_size; j < ny + stencil_size; j++) {
					if (ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1)
						fprintf(f, "%20.10lf\t", B_prev[B_IND(p, i, j, k)]);
				}
				fprintf(f, "\n");
			}
		}
	}
	fclose(f);
	return 0;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S, -t must be declared\n\n");
}

int main(int argc, char **argv)
{
	int i, time_steps;
	in I;
	if (set_parameters(I)) goto error;
	if (read_asc_and_declare_variables(I)) goto error;
	if (do_interpolation(I)) goto error;
	set_arrays(I);
	make_boundary(I);
	I.system_dimension = I.n_cells_multipl * I.nz * I.num_parameters;
	I.system_dimension_with_boundary = I.n_boundary_cells * (I.nz + 2 * I.stencil_size) * I.num_parameters;
	if ((I.B_prev = (double *) malloc(I.system_dimension_with_boundary * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I.B_prev, 0, I.system_dimension_with_boundary * sizeof(double));
	SET_CONDITION(initial, phase_fraction, ascii_map);
	SET_CONDITION(initial, pressure, fixed_value_with_hydrostatic_pressure);
	SET_CONDITION(initial, velocity, mass_is_mooving);
	//SET_CONDITION(initial, velocity, fixed_value);
	set_arrays();
	time_steps = I.end_time / I.dt;
	I.flag_first_time_step = 1;
	for (i = 0; i <= I.time_steps; i++) {
		SET_CONDITION(boundary, velocity, zero_gradient_on_up_and_sides_no_slip_on_low);
		SET_CONDITION(boundary, phase_fraction, fixed_value_on_all);
		SET_CONDITION(boundary, pressure, fixed_value_on_up_and_sides_zero_gradient_on_low);
		if (print_vtk(i) == 1) {
			printf("Error printing vtk file\n");
			goto error;
		} else {
			printf("Result printed to vtk file\n");
		}
		if (create_Ab() == 1) goto error;
		if (i == 0)
			flag_first_time_step = 0;
		if (solve_matrix() == 1) goto error;
		write_B_to_B_prev();
		if (check_conservation_of_mass()) {
			printf("Mass conservation equation failed\n");
		}
		if (i == time_steps) {
			if (print_vtk(i + 1) == 1) {
				printf("Error printing vtk file\n");
				goto error;
			} else {
				printf("Result printed to vtk file\n");
			}
		}
	}
	if (free_massives() == 1) goto error;

	printf("Calculations finished successfully\n");
	return 0;
error:
	printf("Error\n");
	display_usage();
	return 1;
}
