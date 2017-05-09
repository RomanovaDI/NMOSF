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
*/

/*
list of initial conditions
	fixed_value
	ascii_map
	fixed_value_with_hydrostatic_pressure
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include "slu_ddefs.h"

char *map_name;
char *region_map_name;
double hight;
int kx, ky, kz;
int nx, ny, nz;
int n_bl_x, n_bl_y, n_bl_z;
int ncols;
int nrows;
double cellsize;
double *mass;
int *ind;
int *bl_cond;
int n_points, n_cells;
int n_points_multipl, n_cells_multipl;
int *ind_multipl, *ind_cell_multipl;
double interpolation, interpolation_poli;
double *mass1;
double nodata_value;
int *snow_region;
double depth;
double density_snow, density_air;
double pressure_atmosphere;
double g[3];
double k_viscosity_snow, k_viscosity_air;
double *normal, *volume, *area;
double dt, dx, dy, dz;
double *B, *Aelem_csr, *B_prev;
int *Ajptr_csr, *Aiptr_csr;
int num_parameters;
double shear_rate_0,limiting_viscosity_snow, flow_index, yield_stress;
int system_dimension, system_dimension_with_boundary;
int non_zero_elem;
int A_ind_current;
struct timeval tv1,tv2,dtv;
struct timezone tz;
int flag_first_time_step;
double end_time;
int stencil_size;
int *ind_boundary_cells, n_boundary_cells;

int read_asc_and_declare_variables(void);
int do_interpolation(void);
int print_vtk(int n);
int free_massives(void);
double density(int i, int j, int k);
void print_mass(void);
int set_arrays(void);
int conformity(int i, int j);
double velocity_on_face(int p, int i, int j, int k, int s);
double density_on_face(int i, int j, int k, int s);
double strain_rate_on_face(int i, int j, int k, int s, int m, int n);
double shear_rate_on_face(int i, int j, int k, int s);
double phase_fraction_on_face(int i, int j, int k, int s);
double effective_viscosity_on_face(int i, int j, int k, int s);
double pressure_on_face(int i, int j, int k, int s);
double strain_rate(int i, int j, int k, int m, int n);
double shear_rate(int i, int j, int k);
double effective_viscosity(int i, int j, int k);
int DDT_density_velocity_first_combined_VOF(int p, int i, int j, int k);
int DDT_density_velocity_second_combined_VOF(int p, int i, int j, int k);
int DDT_density_velocity_second_separated_VOF(int p, int i, int j, int k);
int DDT_density_velocity_second_separated_FDM(int p, int i, int j, int k);
int DDT_density_velocity_second_combined_FDM(int p, int i, int j, int k);
int DIV_density_velocity_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_forward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_backward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_velocity_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int GRAD_pressure_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int VECT_gravity_force_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int DIV_shear_stress_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int DIV_shear_stress_half_forward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_shear_stress_half_backward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int DIV_grad_pressure_half_forward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_grad_pressure_half_backward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_FDM(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_FDM(int p, int i, int j, int k);
int DDT_density_snow_volume_fraction_first_combined_VOF(int p, int i, int j, int k);
int DDT_density_snow_volume_fraction_second_separated_VOF(int p, int i, int j, int k);
int DDT_density_snow_volume_fraction_second_separated_FDM(int p, int i, int j, int k);
int DDT_density_snow_volume_fraction_second_combined_FDM(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_forward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_backward_euler_second_separated_FDM(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k);
int create_Ab(void);
void display_usage(void);
int A_IND_S(int p, int i, int j, int k, int s);
int A_IND_S_SWITCH(int i, int j, int k, int s);
int A_IND(int p, int i, int j, int k);
int B_IND(int p, int i, int j, int k);
int AREA_IND(int i, int j, int k, int s);
int NORMAL_IND(int p, int i, int j, int k, int s);
int VOLUME_IND(int i, int j, int k);
void print_A_csr(void);
int check_for_corrupt_cell(int i, int j, int k);
int solve_matrix(void);
int write_to_A_csr(int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value);
void time_start(void);
long time_stop(void);
int make_boundary(void);
int SET_initial_CONDITION_velocity_fixed_value(void);
int SET_initial_CONDITION_phase_fraction_ascii_map(void);
int SET_initial_CONDITION_pressure_fixed_value_with_hydrostatic_pressure(void);
int SET_boundary_CONDITION_velocity_zero_gradient_on_up_and_sides_no_slip_on_low(void);
int SET_boundary_CONDITION_phase_fraction_zero_gradient_on_all(void);
int SET_boundary_CONDITION_pressure_fixed_value_on_up_and_sides_zero_gradient_on_low(void);
int write_B_to_B_prev(void);
int write_B_prev_to_B(void);
double velocity(int p, int i, int j, int k);
double phase_fraction(int i, int j, int k);
double pressure(int i, int j, int k);

/* modes of matrix creating */
#define PATTERN 0
#define MATRIX 1

#define DDT(p, i, j, k, object, approximation_order, solution_mode, method) DDT_##object##_##approximation_order##_##solution_mode##_##method(p, i, j, k)
#define DIV(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) DIV_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(p, i, j, k)
#define GRAD(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) GRAD_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(p, i, j, k)
#define VECT(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) VECT_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(p, i, j, k)
#define SET_CONDITION(type, object, mode) SET_##type##_CONDITION_##object##_##mode()

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

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

int A_IND(int p, int i, int j, int k)
{
	return (num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p);
}

int B_IND(int p, int i, int j, int k)
{
	return (num_parameters * ((k + stencil_size) * n_boundary_cells + ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size]) + p);
}

int AREA_IND(int i, int j, int k, int s)
{
	return (6 * ind_cell_multipl[i * ny + j] + s);
}

int NORMAL_IND(int p, int i, int j, int k, int s)
{
	return (18 * ind_cell_multipl[i * ny + j] + s * 3 + p);
}

int VOLUME_IND(int i, int j, int k)
{
	return (ind_cell_multipl[i * ny + j]);
}

double velocity(int p, int i, int j, int k)
{
	return B_prev[B_IND(p, i, j, k)];
}

double phase_fraction(int i, int j, int k)
{
	return B_prev[B_IND(3, i, j, k)];
}

double pressure(int i, int j, int k)
{
	return B_prev[B_IND(4, i, j, k)];
}

int read_asc_and_declare_variables(void)
{
	printf("Reading maps.\n");
	FILE *f = fopen(map_name,"r");
	if (f == NULL) {
		printf("No such file of map\n");
		return 1;
	}

	FILE *f1 = fopen("map.txt", "w");
	int i, j;
	while ((i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);

	f = fopen(region_map_name,"r");
	if (f == NULL) {
		printf("No such file of region\n");
		return 1;
	}

	f1 = fopen("regions_map.txt", "w");
	while ((i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);

	f = fopen("map.txt", "r");
	f1 = fopen("regions_map.txt", "r");
	int err;
	char str[20];

	if ((err = fscanf(f, "%s %d", str, &ncols)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "ncols") != 0) {
		printf("Error file of map: no \"ncols\" tag\n");
		return 1;
	}

	if ((err = fscanf(f, "%s %d", str, &nrows)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "nrows") != 0) {
		printf("Error file of map: no \"nrows\" tag\n");
		return 1;
	}

	double xllcorner;
	if ((err = fscanf(f, "%s %lf", str, &xllcorner)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "xllcorner") != 0) {
		printf("Error file of map: no \"xllcorner\" tag\n");
		return 1;
	}
	
	double yllcorner;
	if ((err = fscanf(f, "%s %lf", str, &yllcorner)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "yllcorner") != 0) {
		printf("Error file of map: no \"yllcorner\" tag\n");
		return 1;
	}
	
	if ((err = fscanf(f, "%s %lf", str, &cellsize)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "cellsize") != 0) {
		printf("Error file of map: no \"cellsize\" tag\n");
		return 1;
	}

	if ((err = fscanf(f, "%s %lf", str, &nodata_value)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "NODATA_value") != 0) {
		printf("Error file of map: no \"NODATA_value\" tag\n");
		return 1;
	}

	int ncols1;
	if ((err = fscanf(f1, "%s %d", str, &ncols1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "ncols") != 0) {
		printf("Error file of region: no \"ncols\" tag\n");
		return 1;
	}

	int nrows1;
	if ((err = fscanf(f1, "%s %d", str, &nrows1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "nrows") != 0) {
		printf("Error file of region: no \"nrows\" tag\n");
		return 1;
	}

	double xllcorner1;
	if ((err = fscanf(f1, "%s %lf", str, &xllcorner1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "xllcorner") != 0) {
		printf("Error file of region: no \"xllcorner\" tag\n");
		return 1;
	}
	
	double yllcorner1;
	if ((err = fscanf(f1, "%s %lf", str, &yllcorner1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "yllcorner") != 0) {
		printf("Error file of region: no \"yllcorner\" tag\n");
		return 1;
	}
	
	double cellsize1;
	if ((err = fscanf(f1, "%s %lf", str, &cellsize1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "cellsize") != 0) {
		printf("Error file of region: no \"cellsize\" tag\n");
		return 1;
	}

	double nodata_value1;
	if ((err = fscanf(f1, "%s %lf", str, &nodata_value1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "NODATA_value") != 0) {
		printf("Error file of region: no \"NODATA_value\" tag\n");
		return 1;
	}

	if (cellsize != cellsize1) {
		printf("Cellsize in both maps need to be the same\n");
		return 1;
	}

	if ((cellsize - (int) cellsize != 0) || (cellsize1 - (int) cellsize1 != 0)) {
		printf("In this vercion the value cellsize need to be integer\n");
		return 1;
	}

	if ((((int) fabs(xllcorner - xllcorner1)) % (int) cellsize) ||
		((int) fabs(yllcorner - yllcorner1) % (int) cellsize) ||
		(fabs(xllcorner - xllcorner1) - (int) fabs(xllcorner - xllcorner1) != 0) ||
		(fabs(yllcorner - yllcorner1) - (int) fabs(yllcorner - yllcorner1) != 0)) {
			printf("Difference between xllcorners of maps and yllcorners need to aligned to cellsize\n");
			return 1;
	}

	if ((mass = (double *) malloc(ncols * nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) mass, 0.0, ncols * nrows * sizeof(double));
	if ((ind = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind, 0, ncols * nrows * sizeof(int));
	if ((bl_cond = (int *) malloc((ncols - 1) * (nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) bl_cond, 0, (ncols - 1) * (nrows - 1) * sizeof(int));

	for (i = 0; i < ncols * nrows; i++) {
		if (err = fscanf(f, "%lf", &mass[i]) == EOF) {
			printf("Error file of map\n");
			return 1;
		}
		ind[i] = -1;
		if (i < (ncols - 1) * (nrows - 1))
			bl_cond[i] = -1;
	}
	fclose(f);
	remove("map.txt");
	
	double *mass_tmp;
	if ((mass_tmp = (double *) malloc(ncols1 * nrows1 * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) mass_tmp, 0.0, ncols1 * nrows1 * sizeof(double));
	for (i = 0; i < ncols1 * nrows1; i++) {
		if (err = fscanf(f1, "%lf", &mass_tmp[i]) == EOF) {
			printf("Error file of region\n");
			return 1;
		}
	}
	remove("regions_map.txt");

	if ((snow_region = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) snow_region, 0, ncols * nrows * sizeof(int));
	double xlucorner = xllcorner - cellsize * nrows;
	double ylucorner = yllcorner;
	double xlucorner1 = xllcorner1 - cellsize1 * nrows1;
	double ylucorner1 = yllcorner1;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if ((xlucorner1 - xlucorner >= 0) && (ylucorner -ylucorner1 >= 0) &&
				(j * cellsize + ylucorner >= ylucorner1) &&
				(i * cellsize + xlucorner >= xlucorner1) &&
				(j * cellsize + ylucorner <= ylucorner1 + ncols1 * cellsize) &&
				(i * cellsize + xlucorner <= xlucorner1 + nrows1 * cellsize)) {
					if (mass_tmp[(i - ((int) (xlucorner1 - xlucorner)) /\
						(int) cellsize) * ncols1 + j -\
						((int) (ylucorner1 - ylucorner)) /\
						(int) cellsize] == 1)
						snow_region[i * ncols + j] = 1;
					else if (mass_tmp[(i - ((int) (xlucorner1 - xlucorner)) /\
						(int) cellsize) * ncols1 + j -\
						((int) (ylucorner1 - ylucorner)) / (int) cellsize] == 0)
							snow_region[i * ncols + j] = 0;
					else
						snow_region[i * ncols + j] = -1;
			} else {
					snow_region[i * ncols + j] = -1;
			}
		}
	}
	free(mass_tmp);

	n_points = 0;
	n_cells = 0;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if ((i - 1 >= 0) && (j + 1 < ncols) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[(i - 1) * ncols + j] != nodata_value) &&
				(mass[i * ncols + j + 1] != nodata_value) &&
				(mass[(i - 1) * ncols + j + 1] != nodata_value)) {
					ind[i * ncols +j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((j + 1 < ncols) && (i + 1 < nrows) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[i * ncols + j + 1] != nodata_value) &&
				(mass[(i + 1) * ncols + j] != nodata_value) &&
				(mass[(i + 1) * ncols + j + 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((i + 1 < nrows) && (j - 1 >= 0) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[(i + 1) * ncols + j] != nodata_value) &&
				(mass[i * ncols + j - 1] != nodata_value) &&
				(mass[(i + 1) * ncols + j - 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[i * ncols + j - 1] != nodata_value) &&
				(mass[(i - 1) * ncols + j] != nodata_value) &&
				(mass[(i - 1) * ncols + j - 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			end_points:
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(ind[i * ncols + j] != -1) &&
				(ind[i * ncols + j - 1] != -1) &&
				(ind[(i - 1) * ncols + j - 1] != -1) &&
				(ind[(i - 1) * ncols + j] != -1)) {
					bl_cond[(i - 1) * (ncols - 1) + j - 1] = n_cells;
					n_cells++;
			}
			if (ind[i * ncols + j] == -1)
				mass[i * ncols + j] = nodata_value;
		}
	}
	nx = (nrows - 1) * kx;
	ny = (ncols - 1) * ky;
	nz = (int) hight / cellsize * kz;
	dx = dy = dz = cellsize;
	printf("Maps entered correctly.\n");
	return 0;
}

int do_interpolation(void)
{
	//printf("Interpolating mesh.\n");
	int i, j, k, l = 0;
	double a, b, d;
	double *c;
	double *e;
	double *f;
	
	n_points_multipl = 0;
	for (i = 0; i < nrows - 1; i++) {
		for (j = 0; j < ncols - 1; j++) {
			if (bl_cond[i * (ncols - 1) + j] != -1) {
				if ((j - 1 < 0) || (bl_cond[i * (ncols - 1) + j - 1] == -1)) {
					n_points_multipl += ky + 1;
					if ((i - 1 >= 0) && (bl_cond[(i - 1) * (ncols - 1) + j] != -1)) {
						n_points_multipl -= 1;
					}
				}
				if ((i - 1 < 0) || (bl_cond[(i - 1) * (ncols - 1) + j] == -1)) {
					n_points_multipl += kx;
					if ((i - 1 >= 0) && (j + 1 < ncols - 1) && (bl_cond[(i - 1) * (ncols - 1) + j + 1] != -1)) {
						n_points_multipl -= 1;
					}
				}
				n_points_multipl += kx * ky;
			}
		}
	}

	if ((ind_multipl = (int *) malloc(((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind_multipl, 0, ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(int));

	int count_ind_multipl = 0;
	if ((mass1 = (double *) malloc(((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) mass1, 0, ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(double));
	for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
		mass1[i] = nodata_value;
	for (i = 0; i < nrows; i++)
		for (j = 0; j < ncols; j++)
			mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky] = mass[i * ncols + j];

	//printf("mass1 after assignment mass\n");
	//for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
	//	printf("%f\n", mass1[i]);

	if ((c = (double *) malloc(ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (double *) malloc(ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (double *) malloc(ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	memset((void *) c, 0, ncols * sizeof(double));
	memset((void *) e, 0, ncols * sizeof(double));
	memset((void *) f, 0, ncols * sizeof(double));

	int ind_start, ind_finish, flag;
	for (i = 0; i < nrows; i++) {
		flag = 0;
		ind_start = 0;
		ind_finish = ncols;
		for (j = 0; j < ncols; j++) {
			if ((ind[i * ncols + j] != -1) && (flag == 0)) {
				ind_start = j;
				flag = 1;
			}
			if ((ind[i * ncols + j] == -1) && (flag == 1)) {
				ind_finish = j;
				flag = 2;
			}
		}
		e[ind_start + 2] = f[ind_start + 2] = 0;
		for (j = ind_start + 3; j < ind_finish - 1; j++) {
			e[j] = - 1. / (4 * e[j - 1] + 1);
			f[j] = (6 * (mass[i * ncols + j + 1] - 2 * mass[i * ncols + j] + mass[i * ncols + j - 1]) / (cellsize * cellsize) - 4 * f[j - 1]) / (4 * e[j - 1] + 1);
		}
		j = ind_finish - 2;
		c[j] = (6 * (mass[i * ncols + j + 1] - 2 * mass[i * ncols + j] + mass[i * ncols + j - 1]) / (cellsize * cellsize) - 4 * f[j]) / (1 + 4 * e[j]);
		for (j = ind_finish - 3; j > ind_start + 1; j--) {
			c[j] = e[j + 1] * c[j + 1] + f[j + 1];
		}
		c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
		for (j = ind_start + 1; j < ind_finish; j++) {
			a = mass[i * ncols + j];
			d = (c[j] - c[j - 1]) / cellsize;
			b = (mass[i * ncols + j] - mass[i * ncols + j - 1]) / cellsize + cellsize * (2 * c[j] + c[j - 1]) / 6;
			if ((ind[i * ncols + j - 1] != -1) && (ind[i * ncols + j] != -1)) {
				for (k = 0; k < ky; k++) {
					mass1[i * ((ncols - 1) * ky + 1) * kx + (j - 1) * ky + k] = a + b * ((double) k * cellsize / (double) ky - cellsize) +
						c[j] * pow((double) k * cellsize / (double) ky - cellsize, 2) / 2 +
						d * pow((double) k * cellsize / (double) ky - cellsize, 3) / 6;
				}
			}
		}
		memset((void *) c, 0, ncols * sizeof(double));
		memset((void *) e, 0, ncols * sizeof(double));
		memset((void *) f, 0, ncols * sizeof(double));
	}

	printf("Interpolation along the y axis have been done\n");
	
	if ((c = (double *) realloc(c, nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (double *) realloc(e, nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (double *) realloc(f, nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	memset((void *) c, 0, nrows * sizeof(double));
	memset((void *) e, 0, nrows * sizeof(double));
	memset((void *) f, 0, nrows * sizeof(double));

	for (j = 0; j < ncols; j++) {
		for (l = 0; l < ky; l++) {
			flag = 0;
			ind_start = 0;
			ind_finish = nrows;
			for (i = 0; i < nrows; i++) {
				if (((ind[i * ncols + j] != -1) && (flag == 0) && (l == 0)) ||
					(((ind[i * ncols + j] != -1) && ((j + 1 == ncols) || (ind[i * ncols + j + 1] != -1))) &&
					 	(flag == 0) && (l != 0))) {
							ind_start = i;
							flag = 1;
				}
				if (((ind[i * ncols + j] == -1) && (flag == 1) && (l == 0)) ||
					(((ind[i * ncols + j] == -1) || ((j + 1 == ncols) || (ind[i * ncols + j + 1] == -1))) &&
					 	(flag == 1) && (l != 0))) {
							ind_finish = i;
							flag = 2;
				}
			}
			e[ind_start + 2] = f[ind_start + 2] = 0;
			for (i = ind_start + 3; i < ind_finish - 1; i++) {
				e[i] = - 1. / (4 * e[i - 1] + 1);
				f[i] = (6 * (mass1[(i + 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l] -
					2 * mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l] +
					mass1[(i - 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l]) /
					(cellsize * cellsize) - 4 * f[i - 1]) / (4 * e[i - 1] + 1);
			}
			i = ind_finish - 2;
			c[i] = (6 * (mass1[(i + 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l] -
				2 * mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l] +
				mass1[(i - 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l]) /
				(cellsize * cellsize) - 4 * f[i]) / (1 + 4 * e[i]);
			for (i = ind_finish - 3; i > ind_start + 1; i--) {
				c[i] = e[i + 1] * c[i + 1] + f[i + 1];
			}
			c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
			for (i = ind_start + 1; i < ind_finish; i++) {
				a = mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l];
				d = (c[i] - c[i - 1]) / cellsize;
				b = (mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l] -
					mass1[(i - 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l]) /
					cellsize + cellsize * (2 * c[i] + c[i - 1]) / 6;
				//printf("j = %d, l = %d, i = %d, a = %f, b = %f, c = %f, d = %f\n", j, l, i, a, b, c[i], d);
				if ((mass1[(i - 1) * ((ncols - 1) * ky + 1) * kx + j * ky + l] != nodata_value) &&
					(mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l] != nodata_value)) {
						for (k = 0; k < kx; k++) {
							mass1[(i - 1) * ((ncols - 1) * ky + 1) * kx + k * ((ncols - 1) * ky + 1) + j * ky + l] =
								a + b * ((double) k * cellsize / (double) kx - cellsize) +
								c[i] * pow((double) k * cellsize / (double) kx - cellsize, 2) / 2 +
								d * pow((double) k * cellsize / (double) kx - cellsize, 3) / 6;
						}
				}
			}
			memset((void *) c, 0, nrows * sizeof(double));
			memset((void *) e, 0, nrows * sizeof(double));
			memset((void *) f, 0, nrows * sizeof(double));
			if (j == ncols - 1) break;
		}
	}

	free(c);
	free(e);
	free(f);

	printf("Interpolation along the x axes have been done\n");

	count_ind_multipl = 0;
	for (i = 0; i < (nrows - 1) * kx + 1; i++) {
		for (j = 0; j < (ncols - 1) * ky + 1; j++) {
			if (mass1[i * ((ncols - 1) * ky + 1) + j] != nodata_value) {
				ind_multipl[i * ((ncols - 1) * ky + 1) + j] = count_ind_multipl;
				count_ind_multipl++;
			} else {
				ind_multipl[i * ((ncols - 1) * ky + 1) + j] = -1;
			}
		}
	}

	if ((ind_cell_multipl = (int *) malloc((nrows - 1) * kx * (ncols - 1) * ky * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind_cell_multipl, 0, (nrows - 1) * kx * (ncols - 1) * ky * sizeof(int));

	for (i = 0; i < (nrows - 1) * kx * (ncols - 1) * ky; i++) {
		ind_cell_multipl[i] = -1;
	}

	for (i = 1; i < (nrows - 1) * kx + 1; i++) {
		for (j = 1; j < (ncols - 1) * ky + 1; j++) {
			if ((ind_multipl[i * ((ncols - 1) * ky + 1) + j] != -1) &&
				(ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j] != -1) &&
				(ind_multipl[i * ((ncols - 1) * ky + 1) + j - 1] != -1) &&
				(ind_multipl[(i - 1) * ((ncols - 1) * ky + 1) + j - 1] != -1)) {
					ind_cell_multipl[(i - 1) * (ncols - 1) * ky + j - 1] = n_cells_multipl;
					n_cells_multipl++;
			}
		}
	}

	dx /= (double) kx;
	dy /= (double) ky;
	dz /= (double) kz;
	return 0;
}

int make_boundary(void)
{
	printf("Make boundary function\n");
	if ((ind_boundary_cells = (int *) malloc((nx + 2 * stencil_size) * (ny + 2 * stencil_size) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind_boundary_cells, -1, (nx + 2 * stencil_size) * (ny + 2 * stencil_size) * sizeof(int));
	int i, j, k, l;
	for (i = stencil_size; i < stencil_size + nx; i++) {
		for (j = stencil_size; j < stencil_size + ny; j++) {
			if (ind_cell_multipl[(i - stencil_size) * ny + j - stencil_size] != -1) {
				for (k = -stencil_size; k <= stencil_size; k++) {
					for (l = -stencil_size; l<= stencil_size; l++) {
						ind_boundary_cells[(i + k) * (ny + 2 * stencil_size) + j + l] = 1;
					}
				}
			}
		}
	}
	k = 0;
	for (i = 0; i < nx + 2 * stencil_size; i++) {
		for (j = 0; j < ny + 2 * stencil_size; j++) {
			if (ind_boundary_cells[i * (ny + 2 * stencil_size) + j] != -1) {
				ind_boundary_cells[i * (ny + 2 * stencil_size) + j] = k;
				k++;
			}
		}
	}
	n_boundary_cells = k;
	FILE *f = fopen("boundary.txt", "w");
	for (i = 0; i < nx + 2 * stencil_size; i++) {
		for (j = 0; j < ny + 2 * stencil_size; j++) {
			fprintf(f, "%d ", ind_boundary_cells[i * (ny + 2 * stencil_size) + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}

int SET_initial_CONDITION_velocity_fixed_value(void)
{
	printf("Set the initial condition for velocity with fixed value for all calculation domain\n");
	int i, j, k, p;
	double fixed_value = 0;
	for (p = 0; p < 3; p++) {
		for (k = 0; k < nz; k++) {
			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {
					if (ind_cell_multipl[i * ny + j] != -1) {
						B_prev[B_IND(p, i, j, k)] = fixed_value;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_phase_fraction_ascii_map(void)
{
	printf("Set the initial condition for phase fraction with values token from the ASCII map\n");
	int i, j, k;
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_boundary_cells[i * ny + j] != -1) {
					if (snow_region[(i / kx) * ncols + j / ky] == 1) {
						if ((double) (k + 1) * (cellsize / kz) <= depth) {
							B_prev[B_IND(3, i, j, k)] = 1;
						} else if (((double) (k + 1) * (cellsize / kz) > depth) && ((double) k * (cellsize / kz)) < depth) {
							B_prev[B_IND(3, i, j, k)] = (depth - (double) k * (cellsize / kz)) / (double) (cellsize / kz);
						} else {
							B_prev[B_IND(3, i, j, k)] = 0;
						}
					} else {
						B_prev[B_IND(3, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_pressure_fixed_value_with_hydrostatic_pressure(void)
{
	printf("Set the initial condition for pressure with fixed value atmosphere pressure in respect hydrostatic pressire\n");
	int i, j, k;
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_boundary_cells[i * ny + j] != -1) {
					if (snow_region[(i / kx) * ncols + j / ky] == 1) {
						if ((double) (k + 1) * (cellsize / kz) <= depth) {
							B_prev[B_IND(4, i, j, k)] = pressure_atmosphere + density_snow * g[2] * (depth - (double) (k * (cellsize / kz)));
							//B_prev[B_IND(4, i, j, k)] = pressure_atmosphere;
						} else if (((double) (k + 1) * (cellsize / kz) > depth) && ((double) k * (cellsize / kz)) < depth) {
							B_prev[B_IND(4, i, j, k)] = pressure_atmosphere + density_snow * g[2] * (depth - (double) (k * (cellsize / kz)));
							//B_prev[B_IND(4, i, j, k)] = pressure_atmosphere;
						} else {
							B_prev[B_IND(4, i, j, k)] = pressure_atmosphere;
						}
					} else {
						B_prev[B_IND(4, i, j, k)] = pressure_atmosphere;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_velocity_zero_gradient_on_up_and_sides_no_slip_on_low(void)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	/*for (p = 0; p < 3; p++) {
		for (k = -stencil_size; k < nz + stencil_size; k++) {
			for (i = -stencil_size; i < nx + stencil_size; i++) {
				for (j = -stencil_size; j < ny + stencil_size; j++) {
					if ((ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1) &&
						((!((i >= 0) && (i < nx) &&
						   (j >= 0) && (j < ny) && 
						   (ind_cell_multipl[i * ny + j] != -1))) ||
						 (k < 0) || (k >= nz))) {
							if (k < 0) {
								B_prev[B_IND(p, i, j, k)] = 0;
							} else if (k >= nz) {
								fl_tmp = 0;
								for (l = 0; l < 2 * stencil_size + 1; l++) {
									for (m = 0; m < 2 * stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < ny) &&
											(ind_cell_multipl[(i + search[l]) * ny + j + search[m]] != -1)) {
												B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i + search[l], j + search[m], nz - 1)];
												fl_tmp = 1;
												break;
										}
									}
									if (fl_tmp) break;
								}
							} else {
								fl_tmp = 0;
								for (l = 0; l < 2 * stencil_size + 1; l++) {
									for (m = 0; m < 2 * stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < ny) &&
											(ind_cell_multipl[(i + search[l]) * ny + j + search[m]] != -1)) {
												B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i + search[l], j + search[m], k)];
												fl_tmp = 1;
												break;
										}
									}
									if (fl_tmp) break;
								}
							}
					}
				}
			}
		}
	}
	for (p = 0; p < 3; p++) {
		for (k = -stencil_size; k < nz + stencil_size; k++) {
			for (i = -stencil_size; i < nx + stencil_size; i++) {
				for (j = -stencil_size; j < ny + stencil_size; j++) {
					if ((ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1) &&
						((!((i >= 0) && (i < nx) &&
						   (j >= 0) && (j < ny) && 
						   (ind_cell_multipl[i * ny + j] != -1))) ||
						 (k < 0) || (k >= nz))) {
							if (k < 0)
								B_prev[B_IND(p, i, j, k)] = 0;
					}
				}
			}
		}
	}*/
	for (p = 0; p < 3; p++) {
		for (k = -stencil_size; k < nz + stencil_size; k++) {
			for (i = -stencil_size; i < nx + stencil_size; i++) {
				for (j = -stencil_size; j < ny + stencil_size; j++) {
					if ((ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1) &&
						((!((i >= 0) && (i < nx) &&
						   (j >= 0) && (j < ny) && 
						   (ind_cell_multipl[i * ny + j] != -1))) ||
						 (k < 0) || (k >= nz))) {
								//printf("p = %d, i = %d, j = %d, k = %d\n", p, i, j, k);
								//printf("B_IND(p, i, j, k) = %d\n", B_IND(p, i, j, k));
							//if (k < 0)
								B_prev[B_IND(p, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_phase_fraction_zero_gradient_on_all(void)
{
	printf("Set the boundary condition for volume snow fraction with zero gradient on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	p = 3;
	for (k = -stencil_size; k < nz + stencil_size; k++) {
		for (i = -stencil_size; i < nx + stencil_size; i++) {
			for (j = -stencil_size; j < ny + stencil_size; j++) {
				if ((ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1) &&
					((!((i >= 0) && (i < nx) &&
					   (j >= 0) && (j < ny) && 
					   (ind_cell_multipl[i * ny + j] != -1))) ||
					 (k < 0) || (k >= nz))) {
						/*if (k < 0) {
							fl_tmp = 0;
							for (l = 0; l < 2 * stencil_size + 1; l++) {
								for (m = 0; m < 2 * stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < ny) &&
										(ind_cell_multipl[(i + search[l]) * ny + j + search[m]] != -1)) {
											B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i + search[l], j + search[m], 0)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						} else if (k >= nz) {
							fl_tmp = 0;
							for (l = 0; l < 2 * stencil_size + 1; l++) {
								for (m = 0; m < 2 * stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < ny) &&
										(ind_cell_multipl[(i + search[l]) * ny + j + search[m]] != -1)) {
											B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i + search[l], j + search[m], nz - 1)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						} else {
							fl_tmp = 0;
							for (l = 0; l < 2 * stencil_size + 1; l++) {
								for (m = 0; m < 2 * stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < ny) &&
										(ind_cell_multipl[(i + search[l]) * ny + j + search[m]] != -1)) {
											B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i + search[l], j + search[m], k)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						}*/
					B_prev[B_IND(p, i, j, k)] = 0;
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_pressure_fixed_value_on_up_and_sides_zero_gradient_on_low(void)
{
	printf("Set the boundary condition for pressure with fixed value of atmospheric pressure on upper and sides boundaries, zero gradient condition for lower boundary\n");
	int i, j, k, p;
	p = 4;
	for (k = nz + stencil_size - 1; k >= - stencil_size; k--) {
		for (i = -stencil_size; i < nx + stencil_size; i++) {
			for (j = -stencil_size; j < ny + stencil_size; j++) {
				if ((ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size] != -1) &&
					((!((i >= 0) && (i < nx) &&
					   (j >= 0) && (j < ny) && 
					   (ind_cell_multipl[i * ny + j] != -1))) ||
					 (k < 0) || (k >= nz))) {
						if (k >= 0)
							B_prev[B_IND(p, i, j, k)] = pressure_atmosphere;
						else
							B_prev[B_IND(p, i, j, k)] = B_prev[B_IND(p, i, j, 0)];
				}
			}
		}
	}
	return 0;
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

int free_massives(void)
{
	free(B);
	free(area);
	free(normal);
	free(volume);
	free(ind_cell_multipl);
	free(mass1);
	free(ind_multipl);
	free(snow_region);
	free(bl_cond);
	free(ind);
	free(mass);
	free(Aelem_csr);
	free(Ajptr_csr);
	free(Aiptr_csr);
	free(B_prev);
	free(ind_boundary_cells);
	return 0;
}

double density(int i, int j, int k)
{
	return phase_fraction(i, j, k) * density_snow + (1 - phase_fraction(i, j, k)) * density_air;
}

void print_mass(void)
{
	int i;
	printf("mass\n");
	for (i = 0; i < ncols * nrows; i++)
		printf("%f\n", mass[i]);
	printf("mass1\n");
	printf("kx = %d\nky = %d\nkz = %d\n", kx, ky, kz);
	for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
		printf("%f\n", mass1[i]);
	return;
}

int set_arrays(void)
{
	printf("Set arrays of volume, normales and areas for mesh\n");
	time_start();
	int i, j, k, l, m;
	if ((volume = (double *) malloc(n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) volume, 0, n_cells_multipl * sizeof(double));
	if ((normal = (double *) malloc(6 * 3 * n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) normal, 0, 6 * 3 * n_cells_multipl * sizeof(double));
	if ((area = (double *) malloc(6 * n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) area, 0, 6 * n_cells_multipl * sizeof(double));
	double a[4], b[4], c[4], p, d;
	k = 0;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (ind_cell_multipl[i * ny + j] != -1) {
				a[1] = cellsize / (double) kx;
				a[2] = 0;
				a[3] = mass1[(i + 1) * (ny + 1) + j] - mass1[i * (ny + 1) + j];
				b[1] = 0;
				b[2] = cellsize / (double) ky;
				b[3] = mass1[i * (ny + 1) + j + 1] - mass1[i * (ny + 1) + j];
				c[1] = cellsize / (double) kx;
				c[2] = - cellsize / (double) ky;
				c[3] = mass1[(i + 1) * (ny + 1) + j] - mass1[i * (ny + 1) + j + 1];
				normal[k * 18 + 4 * 3 + 0] = a[2] * b[3] - a[3] * b[2];
				normal[k * 18 + 4 * 3 + 1] = a[3] * b[1] - a[1] * b[3];
				normal[k * 18 + 4 * 3 + 2] = a[1] * b[2] - a[2] * b[1];
				a[0] = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
				b[0] = sqrt(b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
				c[0] = sqrt(c[1] * c[1] + c[2] * c[2] + c[3] * c[3]);
				p = (a[0] + b[0] + c[0]) / 2;
				area[k * 6 + 4] = sqrt(p * (p - a[0]) * (p - b[0]) * (p - c[0]));
				a[1] = cellsize / (double) kx;
				a[2] = 0;
				a[3] = mass1[(i + 1) * (ny + 1) + j + 1] - mass1[i * (ny + 1) + j + 1];
				b[1] = 0;
				b[2] = cellsize / (double) ky;
				b[3] = mass1[(i + 1) * (ny + 1) + j + 1] - mass1[(i + 1) * (ny + 1) + j];
				normal[k * 18 + 4 * 3 + 0] += a[2] * b[3] - a[3] * b[2];
				normal[k * 18 + 4 * 3 + 1] += a[3] * b[1] - a[1] * b[3];
				normal[k * 18 + 4 * 3 + 2] += a[1] * b[2] - a[2] * b[1];
				a[0] = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
				b[0] = sqrt(b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
				p = (a[0] + b[0] + c[0]) / 2;
				area[k * 6 + 4] += sqrt(p * (p - a[0]) * (p - b[0]) * (p - c[0]));
				d = sqrt(pow(normal[k * 18 + 4 * 3 + 0], 2) + pow(normal[k * 18 + 4 * 3 + 1], 2) + pow(normal[k * 18 + 4 * 3 + 2], 2));
				normal[k * 18 + 4 * 3 + 0] /= d;
				normal[k * 18 + 4 * 3 + 1] /= d;
				normal[k * 18 + 4 * 3 + 2] /= d;
				normal[k * 18 + 5 * 3 + 0] = - normal[k * 18 + 4 * 3 + 0];
				normal[k * 18 + 5 * 3 + 1] = - normal[k * 18 + 4 * 3 + 1];
				normal[k * 18 + 5 * 3 + 2] = - normal[k * 18 + 4 * 3 + 2];
				normal[k * 18 + 2 * 3 + 0] = 0;
				normal[k * 18 + 2 * 3 + 1] = 1;
				normal[k * 18 + 2 * 3 + 2] = 0;
				normal[k * 18 + 3 * 3 + 0] = 0;
				normal[k * 18 + 3 * 3 + 1] = - 1;
				normal[k * 18 + 3 * 3 + 2] = 0;
				normal[k * 18 + 0 * 3 + 0] = 1;
				normal[k * 18 + 0 * 3 + 1] = 0;
				normal[k * 18 + 0 * 3 + 2] = 0;
				normal[k * 18 + 1 * 3 + 0] = - 1;
				normal[k * 18 + 1 * 3 + 1] = 0;
				normal[k * 18 + 1 * 3 + 2] = 0;
				area[k * 6 + 2] = cellsize * cellsize / ((double) kx * (double) kz);
				area[k * 6 + 3] = cellsize * cellsize / ((double) kx * (double) kz);
				area[k * 6 + 0] = cellsize * cellsize / ((double) ky * (double) kz);
				area[k * 6 + 1] = cellsize * cellsize / ((double) ky * (double) kz);
				area[k * 6 + 5] = area[k * 6 + 0];
				volume[k] = cellsize * cellsize * cellsize / ((double) kx * (double) ky * (double) kz);
				k++;
			}
		}
	}
	FILE *f = fopen("volume", "w");
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (ind_cell_multipl[i * ny + j] != -1) {
				fprintf(f, "%20.10lf\t", volume[VOLUME_IND(i, j, k)]);
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen("area", "w");
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (ind_cell_multipl[i * ny + j] != -1) {
				fprintf(f, "%20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf\t", area[AREA_IND(i, j, k, 0)], area[AREA_IND(i, j, k, 1)],
					area[AREA_IND(i, j, k, 2)], area[AREA_IND(i, j, k, 3)], area[AREA_IND(i, j, k, 4)], area[AREA_IND(i, j, k, 5)]);
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen("normal", "w");
	for (i = 0; i < n_cells_multipl; i++) {
		for (j = 0; j < 18; j++) {
			if (normal[i * 18 + j] != normal[j])
				fprintf(f, "pizda\n");
		}
	}
	for (i = 0; i < 18; i++)
		fprintf(f, "%20.10lf ", normal[i]);
	fclose(f);
	printf("Time: %ld\n", time_stop());
	return 0;
}

int conformity(int i, int j)
{
	switch (i + j + i * j) {
		case 0:
			return 0;
		case 1:
			return 1;
		case 2:
			return 2;
		case 3:
			return 3;
		case 5:
			return 4;
		case 8:
			return 5;
	}
}

int check_for_corrupt_cell(int i, int j, int k)
{
	if ((i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0) ||
		(ind_cell_multipl[i * ny + j] == -1)) {
			printf("!!!!!!Operating with corrupt cell!!!!!!\n");
			return 1;
	}
	return 0;
}

int write_to_A_csr(int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value)
{
	if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz) && (ind_cell_multipl[i * ny + j] != -1) && (A_IND_S_SWITCH(i, j, k, s))) {
		int l, m, column, row;
		if (!((i_eq >= 0) && (i_eq < nx) && (j_eq >= 0) && (j_eq < ny) && (k_eq >= 0) && (k_eq < nz) && (ind_cell_multipl[i_eq * ny + j_eq] != -1)))
			printf("something wery strange\n");
		if (s == -1) {
			column = A_IND(p, i, j, k);
		} else {
			column = A_IND_S(p, i, j, k, s);
		}
		row = A_IND(p_eq, i_eq, j_eq, k_eq);
		if (flag_first_time_step) {
			for (l = Aiptr_csr[row]; l <= A_ind_current; l++) {
				if (Ajptr_csr[l] == -1) {
					Aelem_csr[l] = value;
					Ajptr_csr[l] = column;
					break;
				} else if (Ajptr_csr[l] == column) {
					Aelem_csr[l] += value;
					break;
				} else if (Ajptr_csr[l] > column) {
					A_ind_current++;
					if (A_ind_current == non_zero_elem) {
						non_zero_elem++;
						if ((Aelem_csr = (double *) realloc(Aelem_csr, non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((Ajptr_csr = (int *) realloc(Ajptr_csr, non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					for (m = A_ind_current; m > l; m--) {
						Aelem_csr[m] = Aelem_csr[m - 1];
						Ajptr_csr[m] = Ajptr_csr[m - 1];
					}
					Aelem_csr[l] = value;
					Ajptr_csr[l] = column;
					break;
				} else if (l == A_ind_current) {
					A_ind_current++;
					if (A_ind_current == non_zero_elem) {
						non_zero_elem++;
						if ((Aelem_csr = (double *) realloc(Aelem_csr, non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((Ajptr_csr = (int *) realloc(Ajptr_csr, non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					Aelem_csr[A_ind_current] = value;
					Ajptr_csr[A_ind_current] = column;
					break;
				}
			}
		} else {
			for (l = Aiptr_csr[row]; l < Aiptr_csr[row + 1]; l++) {
				if (Ajptr_csr[l] == column) {
					Aelem_csr[l] += value;
					break;
				}
			}
		}
	} else {
		B[A_IND(p_eq, i_eq, j_eq, k_eq)] -= value * B_prev[B_IND(p, i, j, k)];
	}
	return 0;
}

int A_IND_S_SWITCH(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return 0;
			else
				return 1;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return 0;
			else
				return 1;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return 0;
			else
				return 1;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return 0;
			else
				return 1;
		case 4:
			if (k + 1 >= nz)
				return 0;
			else
				return 1;
		case 5:
			if (k - 1 < 0)
				return 0;
			else
				return 1;
	}
	return 1;
}

int A_IND_S(int p, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return A_IND(p, i + 1, j, k);
		case 1:
			return A_IND(p, i - 1, j, k);
		case 2:
			return A_IND(p, i, j + 1, k);
		case 3:
			return A_IND(p, i, j - 1, k);
		case 4:
			return A_IND(p, i, j, k + 1);
		case 5:
			return A_IND(p, i, j, k - 1);
	}
}

double velocity_on_face(int p, int i, int j, int k, int s)
{
	int a, b;
	switch (s) {
		case 0:
			return (velocity(p, i, j, k) + velocity(p, i + 1, j, k)) / 2;
		case 1:
			return (velocity(p, i, j, k) + velocity(p, i - 1, j, k)) / 2;
		case 2:
			return (velocity(p, i, j, k) + velocity(p, i, j + 1, k)) / 2;
		case 3:
			return (velocity(p, i, j, k) + velocity(p, i, j - 1, k)) / 2;
		case 4:
			return (velocity(p, i, j, k) + velocity(p, i, j, k + 1)) / 2;
		case 5:
			return (velocity(p, i, j, k) + velocity(p, i, j, k - 1)) / 2;
	}
}

double density_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (density(i, j, k) + density(i + 1, j, k)) / 2;
		case 1:
			return (density(i, j, k) + density(i - 1, j, k)) / 2;
		case 2:
			return (density(i, j, k) + density(i, j + 1, k)) / 2;
		case 3:
			return (density(i, j, k) + density(i, j - 1, k)) / 2;
		case 4:
			return (density(i, j, k) + density(i, j, k + 1)) / 2;
		case 5:
			return (density(i, j, k) + density(i, j, k - 1)) / 2;
	}
}

double strain_rate_on_face(int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity_on_face(0, i + 1, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / (2 * dx);
		case 1:
			return 0.5 * ((velocity_on_face(0, i, j + 1, k, s) - velocity_on_face(0, i, j - 1, k, s)) / (2 * dy) +
					(velocity_on_face(1, i + 1, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / (2 * dx));
		case 2:
			return 0.5 * ((velocity_on_face(0, i, j, k + 1, s) - velocity_on_face(0, i, j, k - 1, s)) / (2 * dz) +
					(velocity_on_face(2, i + 1, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / (2 * dx));
		case 3:
			return (velocity_on_face(1, i, j + 1, k, s) - velocity_on_face(1, i, j - 1, k, s)) / (2 * dy);
		case 5:
			return 0.5 * ((velocity_on_face(1, i, j, k + 1, s) - velocity_on_face(1, i, j, k - 1, s)) / (2 * dz) +
					(velocity_on_face(2, i, j + 1, k, s) - velocity_on_face(2, i, j - 1, k, s)) / (2 * dy));
		case 8:
			return (velocity_on_face(2, i, j, k + 1, s) - velocity_on_face(2, i, j, k - 1, s)) / (2 * dz);
	}
}

double shear_rate_on_face(int i, int j, int k, int s)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate_on_face(i, j, k, s, m, n) * strain_rate_on_face(i, j, k, s, m, n);
		}
	}
	return x;
}

double phase_fraction_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (phase_fraction(i, j, k) + phase_fraction(i + 1, j, k)) / 2;
		case 1:
			return (phase_fraction(i, j, k) + phase_fraction(i - 1, j, k)) / 2;
		case 2:
			return (phase_fraction(i, j, k) + phase_fraction(i, j + 1, k)) / 2;
		case 3:
			return (phase_fraction(i, j, k) + phase_fraction(i, j - 1, k)) / 2;
		case 4:
			return (phase_fraction(i, j, k) + phase_fraction(i, j, k + 1)) / 2;
		case 5:
			return (phase_fraction(i, j, k) + phase_fraction(i, j, k - 1)) / 2;
	}
}

double effective_viscosity_on_face(int i, int j, int k, int s)
{
	double phase_fraction_on_face_tmp = phase_fraction_on_face(i, j, k, s);
	double shear_rate_on_face_tmp = shear_rate_on_face(i, j, k, s);
	double x = k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= shear_rate_0) {
		x += limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (k_viscosity_snow * pow(shear_rate_on_face_tmp, flow_index - 1) + yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

double effective_viscosity(int i, int j, int k)
{
	double phase_fraction_tmp = phase_fraction(i, j, k);
	double shear_rate_tmp = shear_rate(i, j, k);
	double x = k_viscosity_air * (1 - phase_fraction_tmp);
	if (shear_rate_tmp <= shear_rate_0) {
		x += limiting_viscosity_snow * phase_fraction_tmp;
	} else {
		x += (k_viscosity_snow * pow(shear_rate_tmp, flow_index - 1) + yield_stress / shear_rate_tmp) * phase_fraction_tmp;
	}
	return x;
}

double shear_rate(int i, int j, int k)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate(i, j, k, m, n) * strain_rate(i, j, k, m, n);
		}
	}
	return x;
}

double strain_rate(int i, int j, int k, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity(0, i + 1, j, k) - velocity(0, i - 1, j, k)) / (2 * dx);
		case 1:
			return 0.5 * ((velocity(0, i, j + 1, k) - velocity(0, i, j - 1, k)) / (2 * dy) +
					(velocity(1, i + 1, j, k) - velocity(1, i - 1, j, k)) / (2 * dx));
		case 2:
			return 0.5 * ((velocity(0, i, j, k + 1) - velocity(0, i, j, k - 1)) / (2 * dz) +
					(velocity(2, i + 1, j, k) - velocity(2, i - 1, j, k)) / (2 * dx));
		case 3:
			return (velocity(1, i, j + 1, k) - velocity(1, i, j - 1, k)) / (2 * dy);
		case 5:
			return 0.5 * ((velocity(1, i, j, k + 1) - velocity(1, i, j, k - 1)) / (2 * dz) +
					(velocity(2, i, j + 1, k) - velocity(2, i, j - 1, k)) / (2 * dy));
		case 8:
			return (velocity(2, i, j, k + 1) - velocity(2, i, j, k - 1)) / (2 * dz);
	}
}

double pressure_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (pressure(i, j, k) + pressure(i + 1, j, k)) / 2;
		case 1:
			return (pressure(i, j, k) + pressure(i - 1, j, k)) / 2;
		case 2:
			return (pressure(i, j, k) + pressure(i, j + 1, k)) / 2;
		case 3:
			return (pressure(i, j, k) + pressure(i, j - 1, k)) / 2;
		case 4:
			return (pressure(i, j, k) + pressure(i, j, k + 1)) / 2;
		case 5:
			return (pressure(i, j, k) + pressure(i, j, k - 1)) / 2;
	}
}

int DDT_density_velocity_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	A_value = density(i, j, k) / dt;
	WRITE_TO_A(p, i, j, k, -1);
	B[A_IND(p, i, j, k)] += density(i, j, k) * velocity(p, i, j, k) / dt;
	return 0;
}

int DDT_density_velocity_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_velocity, first, combined, VOF)) return 1;
	return 0;
}

int DDT_density_velocity_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_velocity, second, combined, VOF)) return 1;
	return 0;
}

int DDT_density_velocity_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_velocity, second, combined, VOF)) return 1;
	return 0;
}

int DDT_density_velocity_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_velocity, first, combined, VOF)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_forward_euler, first, combined, VOF)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(p, i, j, k, s) *
				velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)]);
		}
	}
	return 0;
}

int DIV_density_velocity_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(p, i, j, k, s) *
				normal[NORMAL_IND(pr, i, j, k, s)] / ( 2 * 2 * volume[VOLUME_IND(i, j, k)]);
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

int DIV_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] -= (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * density_on_face(i, j, k, s) *
				velocity_on_face(p, i, j, k, s) * velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] +
			 velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]);
		WRITE_TO_A(p, i, j, k, -1);
		WRITE_TO_A(p, i, j, k, s);
	}
	return 0;
}

int DIV_density_velocity_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_forward_euler, second, separated, FDM)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_backward_euler, second, separated, FDM)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_half_forward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pr, i_pr, j_pr, k_pr;
	double d_pr;
	for (pr = 0; pr < 3; pr++) {
		if (pr == 0) {
			d_pr = dx;
			i_pr = 1;
			j_pr = k_pr = 0;
		}
		if (pr == 1) {
			d_pr = dy;
			j_pr = 1;
			i_pr = k_pr = 0;
		}
		if (pr == 2) {
			d_pr = dz;
			k_pr = 1;
			j_pr = i_pr = 0;
		}
		B[A_IND(p, i, j, k)] -= 0.5 *
			(density(i + i_pr, j + j_pr, k + k_pr) * velocity(p, i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_pr, j + j_pr, k + k_pr) -
			 density(i - i_pr, j - j_pr, k - k_pr) * velocity(p, i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_pr, j - j_pr, k - k_pr)) / (2 * d_pr);
	}
	return 0;
}

int DIV_density_velocity_velocity_half_backward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	int pr, i_pr, j_pr, k_pr;
	double d_pr;
	for (pr = 0; pr < 3; pr++) {
		if (pr == 0) {
			d_pr = dx;
			i_pr = 1;
			j_pr = k_pr = 0;
		}
		if (pr == 1) {
			d_pr = dy;
			j_pr = 1;
			i_pr = k_pr = 0;
		}
		if (pr == 2) {
			d_pr = dz;
			k_pr = 1;
			j_pr = i_pr = 0;
		}
		A_value = 0.5 * density(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_pr, j + j_pr, k + k_pr) / (2 * d_pr);
		WRITE_TO_A(p, i + i_pr, j + j_pr, k + k_pr, -1);
		A_value = - 0.5 * density(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_pr, j - j_pr, k - k_pr) / (2 * d_pr);
		WRITE_TO_A(p, i - i_pr, j - j_pr, k - k_pr, -1);
	}
	return 0;
}

int DIV_density_velocity_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, crank_nikolson, second, separated, FDM)) return 1;
	return 0;
}

int GRAD_pressure_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (p == 0) {
		B[A_IND(0, i, j, k)] -= (pressure(i, j, k) - pressure(i - 1, j, k)) / dx;
	}
	if (p == 1) {
		B[A_IND(1, i, j, k)] -= (pressure(i, j, k) - pressure(i, j - 1, k)) / dy;
	}
	if (p == 2) {
		B[A_IND(2, i, j, k)] -= (pressure(i, j, k) - pressure(i, j, k - 1)) / dz;
	}
	return 0;
}

int GRAD_pressure_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, half_forward_euler, first, combined, VOF)) return 1;
	if (GRAD(p, i, j, k, pressure, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int GRAD_pressure_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (p == 0) {
		B[A_IND(0, i, j, k)] -= (pressure(i, j, k) - pressure(i - 1, j, k)) / (2 * dx);
	}
	if (p == 1) {
		B[A_IND(1, i, j, k)] -= (pressure(i, j, k) - pressure(i, j - 1, k)) / (2 * dy);
	}
	if (p == 2) {
		B[A_IND(2, i, j, k)] -= (pressure(i, j, k) - pressure(i, j, k - 1)) / (2 * dz);
	}
	return 0;
}

int GRAD_pressure_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	if (p == 0) {
		A_value = 1 / (2 * dx);
		WRITE_TO_A(4, i, j, k, -1);
		A_value = - 1 / (2 * dx);
		WRITE_TO_A(4, i - 1, j, k, -1);
	}
	if (p == 1) {
		A_value = 1 / (2 * dy);
		WRITE_TO_A(4, i, j, k, -1);
		A_value = - 1 / (2 * dy);
		WRITE_TO_A(4, i, j - 1, k, -1);
	}
	if (p == 2) {
		A_value = 1 / (2 * dz);
		WRITE_TO_A(4, i, j, k, -1);
		A_value = - 1 / (2 * dz);
		WRITE_TO_A(4, i, j, k - 1, -1);
	}
	return 0;
}

int GRAD_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, half_forward_euler, second, combined, VOF)) return 1;
	if (GRAD(p, i, j, k, pressure, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int GRAD_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (p == 0) {
		B[A_IND(0, i, j, k)] -= (pressure(i + 1, j, k) - pressure(i - 1, j, k)) / (4 * dx);
	}
	if (p == 1) {
		B[A_IND(1, i, j, k)] -= (pressure(i, j + 1, k) - pressure(i, j - 1, k)) / (4 * dy);
	}
	if (p == 2) {
		B[A_IND(2, i, j, k)] -= (pressure(i, j, k + 1) - pressure(i, j, k - 1)) / (4 * dz);
	}
	return 0;
}

int GRAD_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	if (p == 0) {
		A_value = 1 / (4 * dx);
		WRITE_TO_A(4, i + 1, j, k, -1);
		A_value = - A_value;
		WRITE_TO_A(4, i - 1, j, k, -1);
	}
	if (p == 1) {
		A_value = 1 / (4 * dy);
		WRITE_TO_A(4, i, j + 1, k, -1);
		A_value = - A_value;
		WRITE_TO_A(4, i, j - 1, k, -1);
	}
	if (p == 2) {
		A_value = 1 / (4 * dz);
		WRITE_TO_A(4, i, j, k + 1, -1);
		A_value = - A_value;
		WRITE_TO_A(4, i, j, k - 1, -1);
	}
	return 0;
}

int GRAD_pressure_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (p == 0) {
		B[A_IND(0, i, j, k)] -= (pressure(i + 1, j, k) - pressure(i - 1, j, k)) / (2 * dx);
	}
	if (p == 1) {
		B[A_IND(1, i, j, k)] -= (pressure(i, j + 1, k) - pressure(i, j - 1, k)) / (2 * dy);
	}
	if (p == 2) {
		B[A_IND(2, i, j, k)] -= (pressure(i, j, k + 1) - pressure(i, j, k - 1)) / (2 * dz);
	}
	return 0;
}

int GRAD_pressure_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, crank_nikolson, second, separated, VOF)) return 1;
	return 0;
}

int GRAD_pressure_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int VECT_gravity_force_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	B[A_IND(p, i, j, k)] += density(i, j, k) * g[p];
	return 0;
}

int VECT_gravity_force_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_forward_euler, first, combined, VOF)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int VECT_gravity_force_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	B[A_IND(p, i, j, k)] += density(i, j, k) * g[p] / 2;
	return 0;
}

int VECT_gravity_force_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	double A_value;
	if (check_for_corrupt_cell(i, j, k)) return 1;
	A_value = - density_snow * g[p] / 2 + density_air * g[p] / 2;
	WRITE_TO_A(3, i, j, k, -1);
	B[A_IND(p, i, j, k)] += density_air * g[p] / 2;
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_forward_euler, second, combined, VOF)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int VECT_gravity_force_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	B[A_IND(p, i, j, k)] += density(i, j, k) * g[p] / 2;
	return 0;
}

int VECT_gravity_force_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	double A_value;
	if (check_for_corrupt_cell(i, j, k)) return 1;
	A_value = - (density_snow - density_air) * g[p] / 2;
	WRITE_TO_A(3, i, j, k, -1);
	B[A_IND(p, i, j, k)] += density_air * g[p] / 2;
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	B[A_IND(p, i, j, k)] += density(i, j, k) * g[p];
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, crank_nikolson, second, separated, VOF)) return 1;
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int DIV_shear_stress_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * strain_rate_on_face(i, j, k, s, p, pr) *
				normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_forward_euler, first, combined, VOF)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int DIV_shear_stress_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] += (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * strain_rate_on_face(i, j, k, s, p, pr) *
				normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, i_p = i, j_p = j, k_p = k, i_pr, j_pr, k_pr;
	double d_p, d_pr, A_value;
	if (p == 0) {
		d_p = dx;
		i_p = i - 1;
	}
	if (p == 1) {
		d_p = dy;
		j_p = j - 1;
	}
	if (p == 2) {
		d_p = dz;
		k_p = k - 1;
	}
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			i_pr = i;
			j_pr = j;
			k_pr = k;
			if (pr == 0) {
				d_pr = dx;
				i_pr = i - 1;
			}
			if (pr == 1) {
				d_pr = dy;
				j_pr = j - 1;
			}
			if (pr == 2) {
				d_pr = dz;
				k_pr = k - 1;
			}
			A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)] * (1 / (2 * d_pr));
			WRITE_TO_A(p, i, j, k, -1);
			WRITE_TO_A(p, i, j, k, s);
			A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)] * (1 / (2 * d_pr));
			WRITE_TO_A(p, i_pr, j_pr, k_pr, -1);
			WRITE_TO_A(p, i_pr, j_pr, k_pr, s);
			A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)] * (1 / (2 * d_p));
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
			A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)] * (1 / (2 * d_p));
			WRITE_TO_A(pr, i_p, j_p, k_p, -1);
			WRITE_TO_A(pr, i_p, j_p, k_p, s);
		}
	}
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_shear_stress_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] += (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				strain_rate_on_face(i, j, k, s, p, pr) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, i_p = 0, j_p = 0, k_p = 0, i_pr, j_pr, k_pr;
	double d_p, d_pr, A_value;
	if (p == 0) {
		d_p = dx;
		i_p = 1;
	}
	if (p == 1) {
		d_p = dy;
		j_p = 1;
	}
	if (p == 2) {
		d_p = dz;
		k_p = 1;
	}
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			i_pr = 0;
			j_pr = 0;
			k_pr = 0;
			if (pr == 0) {
				d_pr = dx;
				i_pr = 1;
			}
			if (pr == 1) {
				d_pr = dy;
				j_pr = 1;
			}
			if (pr == 2) {
				d_pr = dz;
				k_pr = 1;
			}
			A_value = - (area[AREA_IND(i, j, k, s)] / (4 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				normal[NORMAL_IND(pr, i, j, k, s)] / (4 * d_pr);
			WRITE_TO_A(p, i + i_pr, j + j_pr, k + k_pr, -1);
			WRITE_TO_A(p, i + i_pr, j + j_pr, k + k_pr, s);
			A_value = - A_value;
			WRITE_TO_A(p, i - i_pr, j - j_pr, k - k_pr, -1);
			WRITE_TO_A(p, i - i_pr, j - j_pr, k - k_pr, s);
			A_value = - (area[AREA_IND(i, j, k, s)] / (4 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				normal[NORMAL_IND(pr, i, j, k, s)] / (4 * d_p);
			WRITE_TO_A(pr, i + i_p, j + j_p, k + k_p, -1);
			WRITE_TO_A(pr, i + i_p, j + j_p, k + k_p, s);
			A_value = - A_value;
			WRITE_TO_A(pr, i - i_p, j - j_p, k - k_p, -1);
			WRITE_TO_A(pr, i - i_p, j - j_p, k - k_p, s);
		}
	}
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_forward_euler, second, separated, FDM)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_backward_euler, second, separated, FDM)) return 1;
	return 0;
}

int DIV_shear_stress_half_forward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, i_p, j_p, k_p, i_pr, j_pr, k_pr;
	double d;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx * dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy * dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz * dz;
			k_pr = 1;
		}
		B[A_IND(p, i, j, k)] += 0.5 * 2 * 0.5 *
			(effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(p, i + i_pr, j + j_pr, k + k_pr) -
			 2 * effective_viscosity(i, j, k) * velocity(p, i, j, k) +
			 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(p, i - i_pr, j - j_pr, k - k_pr)) / d;
		if (p == pr) {
			B[A_IND(p, i, j, k)] += 0.5 * 2 * 0.5 *
				(effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(p, i + i_pr, j + j_pr, k + k_pr) -
				 2 * effective_viscosity(i, j, k) * velocity(p, i, j, k) +
				 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(p, i - i_pr, j - j_pr, k - k_pr)) / d;
		} else {
			d = 1;
			i_pr = j_pr = k_pr = 0;
			if (pr == 0) {
				d *= dx;
				i_pr = 1;
			}
			if (pr == 1) {
				d *= dy;
				j_pr = 1;
			}
			if (pr == 2) {
				d *= dz;
				k_pr = 1;
			}
			i_p = j_p = k_p = 0;
			if (p == 0) {
				d *= dx;
				i_p = 1;
			}
			if (p == 1) {
				d *= dy;
				j_p = 1;
			}
			if (p == 2) {
				d *= dz;
				k_p = 1;
			}
			B[A_IND(p, i, j, k)] += 0.5 * 2 * 0.5 *
				(effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_p + i_pr, j + j_p + j_pr, k + k_p + k_pr) -
				 effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i - i_p + i_pr, j - j_p + j_pr, k - k_p + k_pr) -
				 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i + i_p - i_pr, j + j_p - j_pr, k + k_p - k_pr) +
				 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_p - i_pr, j - j_p - j_pr, k - k_p - k_pr)) / (4 * d);
		}
	}
	return 0;
}

int DIV_shear_stress_half_backward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, i_p, j_p, k_p, i_pr, j_pr, k_pr;
	double d, A_value;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx * dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy * dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz * dz;
			k_pr = 1;
		}
		A_value = - 0.5 * 2 * 0.5 * effective_viscosity(i + i_pr, j + j_pr, k + k_pr) / d;
		WRITE_TO_A(p, i + i_pr, j + j_pr, k + k_pr, -1);
		A_value = 0.5 * 2 * 0.5 * 2 * effective_viscosity(i, j, k) / d;
		WRITE_TO_A(p, i, j, k, -1);
		A_value = - 0.5 * 2 * 0.5 * effective_viscosity(i - i_pr, j - j_pr, k - k_pr) / d;
		WRITE_TO_A(p, i - i_pr, j - j_pr, k - k_pr, -1);
		if (p == pr) {
			A_value = - 0.5 * 2 * 0.5 * effective_viscosity(i + i_pr, j + j_pr, k + k_pr) / d;
			WRITE_TO_A(p, i + i_pr, j + j_pr, k + k_pr, -1);
			A_value = 0.5 * 2 * 0.5 * 2 * effective_viscosity(i, j, k) / d;
			WRITE_TO_A(p, i, j, k, -1);
			A_value = - 0.5 * 2 * 0.5 * effective_viscosity(i - i_pr, j - j_pr, k - k_pr) / d;
			WRITE_TO_A(p, i - i_pr, j - j_pr, k - k_pr, -1);
		} else {
			d = 1;
			i_pr = j_pr = k_pr = 0;
			if (pr == 0) {
				d *= dx;
				i_pr = 1;
			}
			if (pr == 1) {
				d *= dy;
				j_pr = 1;
			}
			if (pr == 2) {
				d *= dz;
				k_pr = 1;
			}
			i_p = j_p = k_p = 0;
			if (p == 0) {
				d *= dx;
				i_p = 1;
			}
			if (p == 1) {
				d *= dy;
				j_p = 1;
			}
			if (p == 2) {
				d *= dz;
				k_p = 1;
			}
			B[A_IND(p, i, j, k)] += 0.5 * 2 * 0.5 *
				(effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_p + i_pr, j + j_p + j_pr, k + k_p + k_pr) -
				 effective_viscosity(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i - i_p + i_pr, j - j_p + j_pr, k - k_p + k_pr) -
				 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i + i_p - i_pr, j + j_p - j_pr, k + k_p - k_pr) +
				 effective_viscosity(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_p - i_pr, j - j_p - j_pr, k - k_p - k_pr)) / (4 * d);
		}
	}
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, crank_nikolson, second, separated, FDM)) return 1;
	return 0;
}

int DIV_grad_pressure_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_forward_euler, first, combined, VOF)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int DIV_grad_pressure_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i - 1, j, k, s)) * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i, j - 1, k, s)) * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i, j, k - 1, s)) * normal[NORMAL_IND(2, i, j, k, s)] / dz;
	}
	return 0;
}

int DIV_grad_pressure_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		WRITE_TO_A(4, i, j, k, -1);
		WRITE_TO_A(4, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		WRITE_TO_A(4, i - 1, j, k, -1);
		WRITE_TO_A(4, i - 1, j, k, s);
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		WRITE_TO_A(4, i, j, k, -1);
		WRITE_TO_A(4, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		WRITE_TO_A(4, i, j - 1, k, -1);
		WRITE_TO_A(4, i, j - 1, k, s);
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
		WRITE_TO_A(4, i, j, k, -1);
		WRITE_TO_A(4, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
		WRITE_TO_A(4, i, j, k - 1, -1);
		WRITE_TO_A(4, i, j, k - 1, s);
	}
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_grad_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(4, i, j, k)] -=  (area[AREA_IND(i, j, k, s)] / 2) * (
				(pressure_on_face(i + 1, j, k, s) - pressure_on_face(i - 1, j, k, s)) * normal[0, i, j, k, s] / (2 * dx) +
				(pressure_on_face(i, j + 1, k, s) - pressure_on_face(i, j - 1, k, s)) * normal[NORMAL_IND(1, i, j, k, s)] / (2 * dy) +
				(pressure_on_face(i, j, k + 1, s) - pressure_on_face(i, j, k - 1, s)) * normal[NORMAL_IND(2, i, j, k, s)] / (2 * dz));
	}
	return 0;
}

int DIV_grad_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (area[AREA_IND(i, j, k, s)] / 4) * normal[NORMAL_IND(0, i, j, k, s)] / (2 * dx);
		WRITE_TO_A(4, i + 1, j, k, -1);
		WRITE_TO_A(4, i + 1, j, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i - 1, j, k, -1);
		WRITE_TO_A(4, i - 1, j, k, s);
		A_value = (area[AREA_IND(i, j, k, s)] / 4) * normal[NORMAL_IND(1, i, j, k, s)] / (2 * dy);
		WRITE_TO_A(4, i, j + 1, k, -1);
		WRITE_TO_A(4, i, j + 1, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j - 1, k, -1);
		WRITE_TO_A(4, i, j - 1, k, s);
		A_value = (area[AREA_IND(i, j, k, s)] / 4) * normal[NORMAL_IND(2, i, j, k, s)] / (2 * dz);
		WRITE_TO_A(4, i, j, k + 1, -1);
		WRITE_TO_A(4, i, j, k + 1, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j, k - 1, -1);
		WRITE_TO_A(4, i, j, k - 1, s);
	}
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_forward_euler, second, separated, FDM)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_backward_euler, second, separated, FDM)) return 1;
	return 0;
}

int DIV_grad_pressure_half_forward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pr, i_pr, j_pr, k_pr;
	double d;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx * dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy * dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz * dz;
			k_pr = 1;
		}
		B[A_IND(4, i, j, k)] -= 0.5 * (pressure(i + i_pr, j + j_pr, k + k_pr) - 2 * pressure(i, j, k) + pressure(i - i_pr, j - j_pr, k - k_pr)) / d;
	}
	return 0;
}

int DIV_grad_pressure_half_backward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value, d;
	int pr, i_pr, j_pr, k_pr;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx * dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy * dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz * dz;
			k_pr = 1;
		}
		A_value = 1. / (2 * d);
		WRITE_TO_A(4, i + i_pr, j + j_pr, k + k_pr, -1);
		A_value = - 1. / d;
		WRITE_TO_A(4, i, j, k, -1);
		A_value = 1. / (2 * d);
		WRITE_TO_A(4, i - i_pr, j - j_pr, k - k_pr, -1);
	}
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, crank_nikolson, second, separated, FDM)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(4, i, j, k)] -= (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		B[A_IND(4, i, j, k)] -= (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j - 1, k, s)) / dy;
		B[A_IND(4, i, j, k)] -= (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j, k - 1, s)) / dz;
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_forward_euler, first, combined, VOF)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j - 1, k, s)) / dy;
		B[A_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j, k - 1, s)) / dz;
	}
	return 0;
}

int DIV_div_density_velocity_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
		WRITE_TO_A(0, i, j, k, -1);
		WRITE_TO_A(0, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
		WRITE_TO_A(0, i - 1, j, k, -1);
		WRITE_TO_A(0, i - 1, j, k, s);
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
		WRITE_TO_A(1, i, j, k, -1);
		WRITE_TO_A(1, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
		WRITE_TO_A(1, i, j - 1, k, -1);
		WRITE_TO_A(1, i, j - 1, k, s);
		A_value = (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
		WRITE_TO_A(2, i, j, k, -1);
		WRITE_TO_A(2, i, j, k, s);
		A_value = - (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
		WRITE_TO_A(2, i, j, k - 1, -1);
		WRITE_TO_A(2, i, j, k - 1, s);
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(4, i, j, k)] -= (area[AREA_IND(i, j, k, s)] / 2) *
				(((density_on_face(i + 1, j, k, s) * velocity_on_face(pr, i + 1, j, k, s) * velocity_on_face(0, i + 1, j, k, s) -
				   density_on_face(i - 1, j, k, s) * velocity_on_face(pr, i - 1, j, k, s) * velocity_on_face(0, i - 1, j, k, s)) / (2 * dx)) +
				 ((density_on_face(i, j + 1, k, s) * velocity_on_face(pr, i, j + 1, k, s) * velocity_on_face(1, i, j + 1, k, s) -
				   density_on_face(i, j - 1, k, s) * velocity_on_face(pr, i, j - 1, k, s) * velocity_on_face(1, i, j - 1, k, s)) / (2 * dy)) +
				 ((density_on_face(i, j, k + 1, s) * velocity_on_face(pr, i, j, k + 1, s) * velocity_on_face(2, i, j, k + 1, s) -
				   density_on_face(i, j, k - 1, s) * velocity_on_face(pr, i, j, k - 1, s) * velocity_on_face(2, i, j, k - 1, s)) / (2 * dz))) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, pp;
	double A_value, d_pr, i_pr, j_pr, k_pr;
	for (s = 0; s < 6; s++) {
		for (pp = 0; pp < 3; pp++) {
			for (pr = 0; pr < 3; pr++) {
				if (pr == 0) {
					d_pr = dx;
					i_pr = 1;
					j_pr = k_pr = 0;
				}
				if (pr == 1) {
					d_pr = dy;
					j_pr = 1;
					i_pr = k_pr = 0;
				}
				if (pr == 2) {
					d_pr = dz;
					k_pr = 1;
					j_pr = i_pr = 0;
				}
				A_value = (area[AREA_IND(i, j, k, s)] / 4) * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pp, i, j, k, s)] / (2 * d_pr);
				WRITE_TO_A(pp, i + i_pr, j + j_pr, k + k_pr, -1);
				WRITE_TO_A(pp, i + i_pr, j + j_pr, k + k_pr, s);
				A_value = - A_value;
				WRITE_TO_A(pp, i - i_pr, j - j_pr, k - k_pr, -1);
				WRITE_TO_A(pp, i - i_pr, j - j_pr, k - k_pr, s);
			}
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(4, i, j, k)] -= area[AREA_IND(i, j, k, s)] *
				(((density_on_face(i + 1, j, k, s) * velocity_on_face(pr, i + 1, j, k, s) * velocity_on_face(0, i + 1, j, k, s) -
				   density_on_face(i - 1, j, k, s) * velocity_on_face(pr, i - 1, j, k, s) * velocity_on_face(0, i - 1, j, k, s)) / (2 * dx)) +
				 ((density_on_face(i, j + 1, k, s) * velocity_on_face(pr, i, j + 1, k, s) * velocity_on_face(1, i, j + 1, k, s) -
				   density_on_face(i, j - 1, k, s) * velocity_on_face(pr, i, j - 1, k, s) * velocity_on_face(1, i, j - 1, k, s)) / (2 * dy)) +
				 ((density_on_face(i, j, k + 1, s) * velocity_on_face(pr, i, j, k + 1, s) * velocity_on_face(2, i, j, k + 1, s) -
				   density_on_face(i, j, k - 1, s) * velocity_on_face(pr, i, j, k - 1, s) * velocity_on_face(2, i, j, k - 1, s)) / (2 * dz))) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pp, pr, i_pp, j_pp, k_pp, i_pr, j_pr, k_pr;
	double d;
	for (pp = 0; pp < 3; pp++) {
		for (pr = 0; pr < 3; pr++) {
			if (pp == pr) {
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d = dx * dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d = dy * dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d = dz * dz;
					k_pp = 1;
				}
				B[A_IND(4, i, j, k)] -= (density(i + i_pp, j + j_pp, k + k_pp) * velocity(pp, i + i_pp, j+ j_pp, k + k_pp) * velocity(pp, i + i_pp, j + j_pp, k + k_pp) -
						2 * density(i, j, k) * velocity(pp, i, j, k) * velocity(pp, i, j, k) +
						density(i - i_pp, j - j_pp, k - k_pp) * velocity(pp, i - i_pp, j - j_pp, k - k_pp) * velocity(pp, i - i_pp, j - j_pp, k - k_pp)) / d;
			} else {
				d = 1;
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d *= dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d *= dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d *= dz;
					k_pp = 1;
				}
				i_pr = j_pr = k_pr = 0;
				if (pr == 0) {
					d *= dx;
					i_pr = 1;
				}
				if (pr == 1) {
					d *= dy;
					j_pr = 1;
				}
				if (pr == 2) {
					d *= dz;
					k_pr = 1;
				}
				B[A_IND(4, i, j, k)] -=
					(density(i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) * velocity(pp, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) *
					 	velocity(pr, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) -
					 density(i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) * velocity(pp, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) *
					 	velocity(pr, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) -
					 density(i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) * velocity(pp, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) *
					 	velocity(pr, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) +
					 density(i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) * velocity(pp, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) *
					 	velocity(pr, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr)) /
					(4 * d);
			}
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_forward_euler, second, combined, FDM)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_backward_euler, second, combined, FDM)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pp, pr, i_pp, j_pp, k_pp, i_pr, j_pr, k_pr;
	double d;
	for (pp = 0; pp < 3; pp++) {
		for (pr = 0; pr < 3; pr++) {
			if (pp == pr) {
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d = dx * dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d = dy * dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d = dz * dz;
					k_pp = 1;
				}
				B[A_IND(4, i, j, k)] -= 0.5 * (density(i + i_pp, j + j_pp, k + k_pp) * velocity(pp, i + i_pp, j+ j_pp, k + k_pp) * velocity(pp, i + i_pp, j + j_pp, k + k_pp) -
						2 * density(i, j, k) * velocity(pp, i, j, k) * velocity(pp, i, j, k) +
						density(i - i_pp, j - j_pp, k - k_pp) * velocity(pp, i - i_pp, j - j_pp, k - k_pp) * velocity(pp, i - i_pp, j - j_pp, k - k_pp)) / d;
			} else {
				d = 1;
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d *= dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d *= dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d *= dz;
					k_pp = 1;
				}
				i_pr = j_pr = k_pr = 0;
				if (pr == 0) {
					d *= dx;
					i_pr = 1;
				}
				if (pr == 1) {
					d *= dy;
					j_pr = 1;
				}
				if (pr == 2) {
					d *= dz;
					k_pr = 1;
				}
				B[A_IND(4, i, j, k)] -= 0.5 *
					(density(i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) * velocity(pp, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) *
					 	velocity(pr, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) -
					 density(i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) * velocity(pp, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) *
					 	velocity(pr, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) -
					 density(i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) * velocity(pp, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) *
					 	velocity(pr, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) +
					 density(i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) * velocity(pp, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) *
					 	velocity(pr, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr)) /
					(4 * d);
			}
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pp, pr, i_pp, j_pp, k_pp, i_pr, j_pr, k_pr;
	double d, A_value;
	for (pp = 0; pp < 3; pp++) {
		for (pr = 0; pr < 3; pr++) {
			if (pp == pr) {
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d = dx * dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d = dy * dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d = dz * dz;
					k_pp = 1;
				}
				A_value = 0.5 * density(i + i_pp, j + j_pp, k + k_pp) * velocity(pp, i + i_pp, j+ j_pp, k + k_pp) / d;
				WRITE_TO_A(pp, i + i_pp, j + j_pp, k + k_pp, -1);
				A_value = - 0.5 * 2 * density(i, j, k) * velocity(pp, i, j, k) / d;
				WRITE_TO_A(pp, i, j, k, -1);
				A_value = 0,5 * density(i - i_pp, j - j_pp, k - k_pp) * velocity(pp, i - i_pp, j - j_pp, k - k_pp) / d;
				WRITE_TO_A(pp, i - i_pp, j - j_pp, k - k_pp, -1);
			} else {
				d = 1;
				i_pp = j_pp = k_pp = 0;
				if (pp == 0) {
					d *= dx;
					i_pp = 1;
				}
				if (pp == 1) {
					d *= dy;
					j_pp = 1;
				}
				if (pp == 2) {
					d *= dz;
					k_pp = 1;
				}
				i_pr = j_pr = k_pr = 0;
				if (pr == 0) {
					d *= dx;
					i_pr = 1;
				}
				if (pr == 1) {
					d *= dy;
					j_pr = 1;
				}
				if (pr == 2) {
					d *= dz;
					k_pr = 1;
				}
				A_value = 0.5 * density(i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) * velocity(pr, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr) / (4 * d);
				WRITE_TO_A(pp, i + i_pp + i_pr, j + j_pp + j_pr, k + k_pp + k_pr, -1);
				A_value = - 0.5 * density(i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) * velocity(pr, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr) / (4 * d);
				WRITE_TO_A(pp, i + i_pp - i_pr, j + j_pp - j_pr, k + k_pp - k_pr, -1);
				A_value = - 0.5 * density(i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) * velocity(pr, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr) / (4 * d);
				WRITE_TO_A(pp, i - i_pp + i_pr, j - j_pp + j_pr, k - k_pp + k_pr, -1);
				A_value = 0.5 * density(i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) * velocity(pr, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr) / (4 * d);
				WRITE_TO_A(pp, i - i_pp - i_pr, j - j_pp - j_pr, k - k_pp - k_pr, -1);
			}
		}
	}
	return 0;
}

int DDT_density_snow_volume_fraction_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	A_value = 1 / dt;
	WRITE_TO_A(3, i, j, k, -1);
	B[A_IND(3, i, j, k)] += phase_fraction(i, j, k) / dt;
	return 0;
}

int DDT_density_snow_volume_fraction_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_snow_volume_fraction, first, combined, VOF)) return 1;
	return 0;
}

int DDT_density_snow_volume_fraction_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_snow_volume_fraction, first, combined, VOF)) return 1;
	return 0;
}

int DDT_density_snow_volume_fraction_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DDT(p, i, j, k, density_snow_volume_fraction, first, combined, VOF)) return 1;
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_crank_nikolson_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_forward_euler, first, combined, VOF)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_backward_euler, first, combined, VOF)) return 1;
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_forward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(3, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * phase_fraction_on_face(i, j, k, s);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_backward_euler_first_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 4);
		WRITE_TO_A(3, i, j, k, -1);
		WRITE_TO_A(3, i, j, k, s);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(3, i, j, k)] -= (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * phase_fraction_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] +
			 velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (area[AREA_IND(i ,j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] +
			 velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]);
		WRITE_TO_A(3, i, j, k, -1);
		WRITE_TO_A(3, i, j, k, s);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_separated_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, crank_nikolson, second, combined, VOF)) return 1;
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_forward_euler, second, separated, FDM)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, half_backward_euler, second, separated, FDM)) return 1;
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_forward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double d;
	int pr, i_pr, j_pr, k_pr;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz;
			k_pr = 1;
		}
		if ((i == 2) && (j == 3) && (k == 0)) {
			printf("%20.10lf\n", 0.5 * (phase_fraction(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_pr, j + j_pr, k + k_pr) -
				phase_fraction(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_pr, j - j_pr, k - k_pr)) / (2 * d));
		}
		B[A_IND(3, i, j, k)] -= 0.5 * (phase_fraction(i + i_pr, j + j_pr, k + k_pr) * velocity(pr, i + i_pr, j + j_pr, k + k_pr) -
				phase_fraction(i - i_pr, j - j_pr, k - k_pr) * velocity(pr, i - i_pr, j - j_pr, k - k_pr)) / (2 * d);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_half_backward_euler_second_separated_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value, d;
	int pr, i_pr, j_pr, k_pr;
	for (pr = 0; pr < 3; pr++) {
		i_pr = j_pr = k_pr = 0;
		if (pr == 0) {
			d = dx;
			i_pr = 1;
		}
		if (pr == 1) {
			d = dy;
			j_pr = 1;
		}
		if (pr == 2) {
			d = dz;
			k_pr = 1;
		}
		A_value = 0.5 * velocity(pr, i + i_pr, j + j_pr, k + k_pr) / (2 * d);
		WRITE_TO_A(3, i + i_pr, j + j_pr, k + k_pr, -1);
		A_value = - 0.5 * velocity(pr, i - i_pr, j - j_pr, k - k_pr) / (2 * d);
		WRITE_TO_A(3, i - i_pr, j - j_pr, k - k_pr, -1);
	}
	return 0;
}

int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_combined_FDM(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, crank_nikolson, second, separated, FDM)) return 1;
	return 0;
}

int write_B_to_B_prev(void)
{
	int i, j, k, p;
	for (p = 0; p < num_parameters; p++) {
		for (k = 0; k < nz; k++) {
			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {
					if (ind_cell_multipl[i * ny + j] != -1)
						B_prev[B_IND(p, i, j, k)] = B[A_IND(p, i, j, k)];
				}
			}
		}
	}
	return 0;
}

int write_B_prev_to_B(void)
{
	int i, j, k, p;
	for (p = 0; p < num_parameters; p++) {
		for (k = 0; k < nz; k++) {
			for (i = 0; i < nx; i++) {
				for (j = 0; j < ny; j++) {
					if (ind_cell_multipl[i * ny + j] != -1)
						B[A_IND(p, i, j, k)] = B_prev[B_IND(p, i, j, k)];
				}
			}
		}
	}
	return 0;
}

int create_Ab(void)
{
	printf("Creating matrix A and vector B function\n");
	time_start();
	int i, j, k, p;
	if (flag_first_time_step) {
		if ((B = (double *) malloc(system_dimension * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) B, 0, system_dimension * sizeof(double));
		int num_el_in_row = 100; // it is inaccurate value
		non_zero_elem = num_el_in_row * system_dimension;
		if ((Aelem_csr = (double *) malloc(non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) Aelem_csr, 0, non_zero_elem * sizeof(double));
		if ((Ajptr_csr = (int *) malloc(non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) Ajptr_csr, -1, non_zero_elem * sizeof(int));
		if ((Aiptr_csr = (int *) malloc((system_dimension + 1) * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) Aiptr_csr, -1, (system_dimension + 1) * sizeof(int));
	} else {
		memset((void *) B, 0, system_dimension * sizeof(double));
		memset((void *) Aelem_csr, 0, non_zero_elem * sizeof(double));
	}

/* creating matrix */
	if (flag_first_time_step)
		A_ind_current = 0;
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1) {
					/* momentum equation */
					for (p = 0; p < 3; p++) {
						if (flag_first_time_step) {
							if (Ajptr_csr[A_ind_current] != -1)
								A_ind_current++;
							Aiptr_csr[A_IND(p, i, j, k)] = A_ind_current;
						}
						if (DDT(p, i, j, k, density_velocity, second, combined, FDM)) return 1;
						if (DIV(p, i, j, k, density_velocity_velocity, crank_nikolson, second, combined, FDM)) return 1;
						if (VECT(p, i, j, k, gravity_force, crank_nikolson, second, combined, FDM)) return 1;
						if (GRAD(p, i, j, k, pressure, crank_nikolson, second, combined, FDM)) return 1;
						if (DIV(p, i, j, k, shear_stress, crank_nikolson, second, combined, FDM)) return 1;
					}
					/* transport equation for snow volume fraction */
					p = 3;
					if (flag_first_time_step) {
						if (Ajptr_csr[A_ind_current] != -1)
							A_ind_current++;
						Aiptr_csr[A_IND(p, i, j, k)] = A_ind_current;
					}
					if (DDT(p, i, j, k, density_snow_volume_fraction, second, combined, FDM)) return 1;
					if (DIV(p, i, j, k, density_snow_volume_fraction_velocity, crank_nikolson, second, combined,FDM)) return 1;
					/* poisson equation for pressure */
					p = 4;
					if (flag_first_time_step) {
						if (Ajptr_csr[A_ind_current] != -1)
							A_ind_current++;
						Aiptr_csr[A_IND(p, i, j, k)] = A_ind_current;
					}
					if (DIV(p, i, j, k, grad_pressure, crank_nikolson, second, combined, FDM)) return 1;
					if (DIV(p, i, j, k, div_density_velocity_velocity, crank_nikolson, second, combined, FDM)) return 1;
				}
			}
		}
	}
	if (flag_first_time_step) {
		non_zero_elem = A_ind_current + 1;
		Aiptr_csr[system_dimension] = non_zero_elem;
		if ((Aelem_csr = (double *) realloc(Aelem_csr, non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		if ((Ajptr_csr = (int *) realloc(Ajptr_csr, non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
	}
	//if (flag_first_time_step)
	//	print_A_csr();
	printf("Time: %ld\n", time_stop());
	return 0;
}

void print_A_csr(void)
{
	printf("Print matrix A in CSR format function\n");
	time_start();
	FILE *f;
	int i, j, k, fl_tmp;
	if ((f = fopen("A_csr.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", system_dimension);
	double max = abs(Aelem_csr[0]);
	for (i = 0; i < non_zero_elem; i++) {
		if (max < abs(Aelem_csr[i]))
			max = abs(Aelem_csr[i]);
		j = 0;
		while (!((i >= Aiptr_csr[j]) && (i < Aiptr_csr[j + 1])))
			j++;
		fprintf(f, "%20.10lf\t%10d\t%10d\n", Aelem_csr[i], j, Ajptr_csr[i]);
	}
	fclose(f);
	if ((f = fopen("A_B.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	for (i = 0; i < system_dimension; i++) {
		for (j = 0; j < system_dimension; j++) {
			fl_tmp = 1;
			for (k = Aiptr_csr[i]; k < Aiptr_csr[i + 1]; k++) {
				if (Ajptr_csr[k] == j) {
					fprintf(f, "%20.10lf\t", Aelem_csr[k]);
					fl_tmp = 0;
					break;
				}
				if (fl_tmp) {
					fprintf(f, "%20.10lf\t", 0.0);
				}
			}
		}
		fprintf(f, "\t\t%20.10lf\n", B[i]);
	}
	fclose(f);
	if ((f = fopen("A_pattern.dat","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", system_dimension);
	for (i = 0; i < non_zero_elem; i++) {
		j = 0;
		while (!((i >= Aiptr_csr[j]) && (i <= Aiptr_csr[j + 1])))
			j++;
		//fprintf(f, "%20.10lf\t%10d\t%10d\n", Aelem_csr[i] / max, system_dimension - 1 - j, Ajptr_csr[i]);
		fprintf(f, "%20.10lf\t%10d\t%10d\n", 1., system_dimension - 1 - j, Ajptr_csr[i]);
	}
	fclose(f);
	int *A_pattern = (int *) malloc(system_dimension * system_dimension * sizeof(int));
	memset(A_pattern, 0, system_dimension * system_dimension * sizeof(int));
	for (i = 0; i < non_zero_elem; i++) {
		j = 0;
		while (!((i >= Aiptr_csr[j]) && (i <= Aiptr_csr[j + 1])))
			j++;
		A_pattern[(system_dimension - 1 - j) * system_dimension + Ajptr_csr[i]] = 1;
	}
	if ((f = fopen("A_pattern_matrix.dat","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", system_dimension);
	for (i = 0; i < system_dimension; i++) {
		for (j = 0; j < system_dimension; j++) {
			fprintf(f, "%d\t", A_pattern[i * system_dimension + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	if ((f = fopen("A_elem_jprt.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", system_dimension);
	for (i = 0; i < non_zero_elem; i++) {
		fprintf(f, "%10d\t%20.10lf\t%10d\n", i, Aelem_csr[i], Ajptr_csr[i]);
	}
	fclose(f);
	if ((f = fopen("A_iptr.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", system_dimension);
	for (i = 0; i < system_dimension; i++) {
		fprintf(f, "%10d\n", Aiptr_csr[i]);
	}
	fclose(f);
	printf("Time: %ld\n", time_stop());
	return;
}

int solve_matrix(void)
{
	printf("Solving matrix\n");
	SuperMatrix A_csr;
	NCformat *Astore;
//	double   *a;
//	int      *asub, *xa;
	int      *perm_c; /* column permutation vector */
	int      *perm_r; /* row permutations from partial pivoting */
	SuperMatrix L;      /* factor L */
	SCformat *Lstore;
	SuperMatrix U;      /* factor U */
	NCformat *Ustore;
	SuperMatrix B_csr;
//	int      nrhs, ldx, info, m, n, nnz;
	int      info, i, j, nrhs;
//	double   *xact, *rhs;
//	double   *xact;
	mem_usage_t   mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Enter main()");
#endif

	/* Set the default input options:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
	 */

	set_default_options(&options);

	/* Read the matrix in Harwell-Boeing format. */
//	dreadhb(&m, &n, &nnz, &a, &asub, &xa);
	dCreate_CompCol_Matrix(&A_csr, system_dimension, system_dimension, non_zero_elem, Aelem_csr, Ajptr_csr, Aiptr_csr, SLU_NR, SLU_D, SLU_GE);
	Astore = A_csr.Store;
	printf("Dimension %dx%d; # nonzeros %d\n", A_csr.nrow, A_csr.ncol, Astore->nnz);
	nrhs   = 1;
//	if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B_csr, system_dimension, nrhs, B, system_dimension, SLU_DN, SLU_D, SLU_GE);
//	xact = doubleMalloc(n * nrhs);
//	ldx = n;
//	dGenXtrue(n, nrhs, xact, ldx);
//	dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
	if ( !(perm_c = intMalloc(system_dimension)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(perm_r = intMalloc(system_dimension)) ) ABORT("Malloc fails for perm_r[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, &A_csr, perm_c, perm_r, &L, &U, &B_csr, &stat, &info);
	if ( info == 0 ) {
		/* This is how you could access the solution matrix. */
		double *sol = (double*) ((DNformat*) B_csr.Store)->nzval; 
		for (i = 0; i < n_cells_multipl * nz * num_parameters; i++) {
			B[i] = sol[i];
		}
		/* Compute the infinity norm of the error. */
//		dinf_norm_error(nrhs, &B, xact);
		dinf_norm_error(nrhs, &B_csr, sol);
		Lstore = (SCformat *) L.Store;
		Ustore = (NCformat *) U.Store;
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - system_dimension);
		printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - system_dimension)/non_zero_elem);
		dQuerySpace(&L, &U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	} else {
		printf("dgssv() error returns INFO= %d\n", info);
		if ( info <= system_dimension ) { /* factorization completes */
			dQuerySpace(&L, &U, &mem_usage);
			printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		}
	}
	if ( options.PrintStat ) StatPrint(&stat);
//	for (i = 0; i < n_cells_multipl * nz; i++) {
//		velocity[i * 3 + 0] = sol[i * num_parameters + 0];
//		velocity[i * 3 + 1] = sol[i * num_parameters + 1];
//		velocity[i * 3 + 3] = sol[i * num_parameters + 2];
//		phase_fraction[i] = sol[i * num_parameters + 3];
//		pressure[i] = sol[i * num_parameters + 4];
//	}
	printf("Matrix solved\n");
	
	StatFree(&stat);
//	SUPERLU_FREE (rhs);
//	SUPERLU_FREE (xact);
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);
//	Destroy_CompCol_Matrix(&A_csr);
//	Destroy_SuperMatrix_Store(&B_csr);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Exit main()");
#endif
	printf("End of function solve matrix\n");
	return 0;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S, -t must be declared\n\n");
	printf("-m\tmap\tASCII file of map\n");
	printf("-r\trerion\tASCII file of region\n");
	printf("-H\thight\thight of calculation area (double)\n");
	printf("-D\tdepth\tdepth of snow cover (double)\n");
	printf("-x\tkx\treduction ratio in x direction (int)\n");
	printf("-y\tky\treduction ratio in y direction (int)\n");
	printf("-z\tkz\treduction ratio in z direction (int)\n");
	printf("-a\tsnow density air density\tdensity of snow and air (double)\n");
	printf("-s\tsnow density air density\tdensity of snow and air (double)\n");
	printf("-p\tpressure\tatmosphere pressure\n");
	//printf("-v\tsnow kinematic viscosity air kinematic viscosity\tkinematic viscisity of snow and air (double)\n");
	//printf("-v\tsnow viscosity air viscosity\tviscisity of snow and air (double)\n");
	printf("-v\tair viscosity\tshear viscosity of the air for Newtonian model (double)\n");
	printf("-k\tsnow consistency index\tconsistency index for the constitutive equation of the Herschel-Bulkley model for snow (double)\n");
	printf("-i\tflow index\tflow index for the constitutive equation of the Herschel-Bulkley model for snow (double)\n");
	printf("-l\tyield stress\tyield stress for the constitutive equation of the Herschel-Bulkley model for snow (double)\n");
	printf("-S\tlimiting shear rate\tlimiting shear rate for the effective viscosity to use Herschel-Bulkley model as a generalized Newtonian fluid model (double)\n");
	printf("-t\ttime interval\tvalue of the time interval for calculations\n");
	printf("-h\tdisplay usage\n");
}

int main(int argc, char **argv)
{
	int opt = 0, i, time_steps;
	g[0] = g[1] = 0;
	g[2] = 9,81;
	//g[2] = 0;
	stencil_size = 2;
	num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	static const char *optString = "m:r:H:D:x:y:z:a:s:p:v:k:i:l:S:t:h?";
	while ((opt = getopt(argc, argv, optString)) != -1) {
		switch (opt) {
			case 'm':
				map_name = optarg;
				break;
			case 'r':
				region_map_name = optarg;
				break;
			case 'H':
				hight = atof(optarg);
				break;
			case 'D':
				depth = atof(optarg);
				break;
			case 'x':
				kx = atoi(optarg);
				break;
			case 'y':
				ky = atoi(optarg);
				break;
			case 'z':
				kz = atoi(optarg);
				break;
			case 'a':
				density_snow = atof(optarg);
				break;
			case 's':
				density_air = atof(optarg);
				break;
			case 'p':
				pressure_atmosphere = atof(optarg);
				break;
			case 'v':
				k_viscosity_air = atof(optarg);
				break;
			case 'k':
				k_viscosity_snow = atof(optarg);
				break;
			case 'i':
				flow_index = atof(optarg);
				break;
			case 'l':
				yield_stress = atof(optarg);
				break;
			case 'S':
				shear_rate_0 = atof(optarg);
				limiting_viscosity_snow = k_viscosity_snow * pow(shear_rate_0, flow_index - 1) + yield_stress / shear_rate_0;
				break;
			case 't':
				end_time = atof(optarg);
				break;
			case 'h':
			case '?':
				display_usage();
				return 1;
		}
	}
	if (argc != 33) {
		printf("Not enouth arguments\n");
		goto error;
	}

	if (read_asc_and_declare_variables() == 1) goto error;
	if (do_interpolation() == 1) goto error;
	set_arrays();
	make_boundary();
	system_dimension = n_cells_multipl * nz * num_parameters;
	system_dimension_with_boundary = n_boundary_cells * (nz + 2 * stencil_size) * num_parameters;
	if ((B_prev = (double *) malloc(system_dimension_with_boundary * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) B_prev, 0, system_dimension_with_boundary * sizeof(double));
	SET_CONDITION(initial, velocity, fixed_value);
	SET_CONDITION(initial, phase_fraction, ascii_map);
	SET_CONDITION(initial, pressure, fixed_value_with_hydrostatic_pressure);
	set_arrays();
	dt = 0.01;//we need to set dt!!!
	flag_first_time_step = 1;
	time_steps = end_time / dt;
	for (i = 0; i <= time_steps; i++) {
		SET_CONDITION(boundary, velocity, zero_gradient_on_up_and_sides_no_slip_on_low);
		SET_CONDITION(boundary, phase_fraction, zero_gradient_on_all);
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
