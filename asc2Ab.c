#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
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
double *phase_fraction;
double *velocity;
double *pressure;
double pressure_atmosphere;
double g[3];
double k_viscosity_snow, k_viscosity_air;
double *normal, *volume, *area;
double dt, dx, dy, dz;
double *A, *B, *Aelem_csr, *Aelem_csc;
int *Ajptr_csr, *Aiptr_csr, *A_stencil, *Aiptr_csc, *Ajptr_csc;
int num_parameters;
double shear_rate_0,limiting_viscosity_snow, flow_index, yield_stress;
int phase_fraction_len, velocity_len, pressure_len, volume_len, normal_len, area_len, A_len, B_len, mass_len, ind_len, bl_cond_len, snow_region_len, ind_multipl_len, mass1_len, ind_cell_multipl_len;
int system_dimension;
int non_zero_elem;

int read_asc_and_declare_variables(void);
int do_interpolation(void);
int print_vtk(int n);
int free_massives(void);
double density(int i, int j, int k);
void print_mass(void);
int set_arrays(void);
int conformity(int i, int j);
void DDT_density_velocity(int i, int j, int k);
double velocity_on_face(int p, int i, int j, int k, int s);
double density_on_face(int i, int j, int k, int s);
void DIV_density_velocity_velocity_half_forward_euler(int i, int j, int k);
int A_IND_S(int p_eq, int i_eq, int j_eq, int k_eq, int p_el, int i_el, int j_el, int k_el, int s);
int A_IND_S_SWITCH(int i, int j, int k, int s);
void DIV_density_velocity_velocity_half_backward_euler(int i, int j, int k);
void DIV_density_velocity_velocity_crank_nikolson(int i, int j, int k);
void GRAD_pressure_crank_nikolson(int i, int j, int k);
void GRAD_pressure_half_forward_euler(int i, int j, int k);
void GRAD_pressure_half_backward_euler(int i, int j, int k);
void VECT_gravity_force_half_forward_euler(int i, int j, int k);
void VECT_gravity_force_half_backward_euler(int i, int j, int k);
void VECT_gravity_force_crank_nikolson(int i, int j, int k);
double strain_rate_on_face(int i, int j, int k, int s, int m, int n);
double shear_rate_on_face(int i, int j, int k, int s);
double phase_fraction_on_face(int i, int j, int k, int s);
double effective_viscosity_on_face(int i, int j, int k, int s);
void DIV_shear_stress_forward_euler(int i, int j, int k);
double pressure_on_face(int i, int j, int k, int s);
void DIV_grad_pressure_crank_nikolson(int i, int j, int k);
void DIV_grad_pressure_half_forward_euler(int i, int j, int k);
void DIV_grad_pressure_half_backward_euler(int i, int j, int k);
void DIV_div_density_velocity_velocity_crank_nikolson(int i, int j, int k);
void DIV_div_density_velocity_velocity_half_forward_euler(int i, int j, int k);
void DIV_div_density_velocity_velocity_half_backward_euler(int i, int j, int k);
void DDT_density_snow_volume_fraction(int i, int j, int k);
void DIV_density_snow_volume_fraction_velocity_crank_nikolson(int i, int j, int k);
void DIV_density_snow_volume_fraction_velocity_half_forward_euler(int i, int j, int k);
void DIV_density_snow_volume_fraction_velocity_half_backward_euler(int i, int j, int k);
int create_Ab(void);
void display_usage(void);
int A_IND(int p_eq, int i_eq, int j_eq, int k_eq, int p_el, int i_el, int j_el, int k_el);
int B_IND(int p, int i, int j, int k);
int PRESS_IND(int i, int j, int k);
int VEL_IND(int p, int i, int j, int k);
int AREA_IND(int i, int j, int k, int s);
int NORMAL_IND(int p, int i, int j, int k, int s);
int VOLUME_IND(int i, int j, int k);
int PHASE_IND(int i, int j, int k);
void print_A_csr(void);
int A_csr_func(void);
int A_csc_func(void);
void printf_Ab(void);
int check_for_corrupt_cell(int i, int j, int k);
int solve_matrix(void);

#define DDT(i, j, k, object) DDT_##object(i, j, k);
#define DIV(i, j, k, object, numerical_scheme) DIV_##object##_##numerical_scheme(i, j, k);
#define BOUNDARY_CONDITION(p, i, j, k, s, object, mode, numerical_scheme) BOUNDARY_CONDITION_##object##_##mode##_##numerical_scheme(p, i, j, k, s)
#define DDX(p, i, j, k, s, mode, numerical_scheme) DDX_##mode##_##numerical_scheme(p, i, j, k, s)
#define GRAD(i, j, k, mode, numerical_scheme) GRAD_##mode##_##numerical_scheme(i, j, k)
#define VECT(i, j, k, mode, numerical_scheme) VECT_##mode##_##numerical_scheme(i, j, k)
#define DDY(p, i, j, k, s, mode, numerical_scheme) DDY_##mode##_##numerical_scheme(p, i, j, k, s)
#define DDZ(p, i, j, k, s, mode, numerical_scheme) DDZ_##mode##_##numerical_scheme(p, i, j, k, s)

int A_IND(int p_eq, int i_eq, int j_eq, int k_eq, int p_el, int i_el, int j_el, int k_el)
{
	int x = ((num_parameters * (k_eq * n_cells_multipl + ind_cell_multipl[i_eq * ny + j_eq]) + p_eq) *
		n_cells_multipl * nz * num_parameters + num_parameters * (k_el * n_cells_multipl + ind_cell_multipl[i_el * ny + j_el]) + p_el);
	if (x >= A_len)
		printf("Error: out of memory\n");
	if ((x == 3) || (x == 2))
		printf("!!!!!!!!!!!!!\n");
	return x;
}

int B_IND(int p, int i, int j, int k)
{
	int x = (num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p);
	if (x >= B_len)
		printf("Error: out of memory\n");
	return x;
}

int PRESS_IND(int i, int j, int k)
{
	int x = (k * n_cells_multipl + ind_cell_multipl[i * ny + j]);
	if (x >= pressure_len)
		printf("Error: out of memory\n");
	return x;
}

int VEL_IND(int p, int i, int j, int k)
{
	int x = (3 * PRESS_IND(i, j, k) + p);
	if (x >= velocity_len)
		printf("Error: out of memory\n");
	return x;
}

int AREA_IND(int i, int j, int k, int s)
{
	int x = (6 * ind_cell_multipl[i * ny + j] + s);
	if (x >= area_len)
		printf("Error: out of memory\n");
	return x;
}

int NORMAL_IND(int p, int i, int j, int k, int s)
{
	int x = (18 * ind_cell_multipl[i * ny + j] + s * 3 + p);
	if (x >= normal_len)
		printf("Error: out of memory\n");
	return x;
}

int VOLUME_IND(int i, int j, int k)
{
	int x = (ind_cell_multipl[i * ny + j]);
	if (x >= volume_len)
		printf("Error: out of memory\n");
	return x;
}

int PHASE_IND(int i, int j, int k)
{
	int x = (k * n_cells_multipl + ind_cell_multipl[i * ny + j]);
	if (x >= phase_fraction_len)
		printf("Error: out of memory\n");
	return x;
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

	mass_len = ncols * nrows;
	if ((mass = (double *) malloc(ncols * nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) mass, 0, ncols * nrows * sizeof(double));
	ind_len = ncols * nrows;
	if ((ind = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind, 0, ncols * nrows * sizeof(int));
	bl_cond_len = (ncols - 1) * (nrows - 1);
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
	memset((void *) mass_tmp, 0, ncols1 * nrows1 * sizeof(double));
	for (i = 0; i < ncols1 * nrows1; i++) {
		if (err = fscanf(f1, "%lf", &mass_tmp[i]) == EOF) {
			printf("Error file of region\n");
			return 1;
		}
	}
	remove("regions_map.txt");

	snow_region_len = ncols * nrows;
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

	ind_multipl_len = ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1);
	if ((ind_multipl = (int *) malloc(((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) ind_multipl, 0, ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(int));

	int count_ind_multipl = 0;
	mass1_len = ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1);
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

	ind_cell_multipl_len = (nrows - 1) * kx * (ncols - 1) * ky;
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

int print_vtk(int n)
{
	int nn;
	int i = 0, j, k, a;
	while (nn > 0) {
		nn /= 10;
		i++;
	}
	char file_name[7 + i];
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
	fprintf(f, "VECTOR velocity double\n");
	for (i = 0; i < n_cells_multipl * nz; i++)
		fprintf(f, "%20.10f\t%20.10f\t%20.10f\n", velocity[i * 3 + 0], velocity[i * 3 + 1], velocity[i * 3 + 2]);
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (i = 0; i < n_cells_multipl * nz; i++)
		fprintf(f, "%20.10f\n", pressure[i]);
	fprintf(f, "SCALARS snow_volume_fraction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (i = 0; i < n_cells_multipl * nz; i++)
		fprintf(f, "%20.10f\n", phase_fraction[i]);
	fclose(f);
	return 0;
}

int free_massives(void)
{
	printf("free(B);\n");
	free(B);
	printf("free(A);\n");
	free(A);
	printf("free(area);\n");
	free(area);
	printf("free(normal);\n");
	free(normal);
	printf("free(volume);\n");
	free(volume);
	printf("free(pressure);\n");
	free(pressure);
	printf("free(velocity);\n");
	free(velocity);
	printf("free(phase_fraction);\n");
	free(phase_fraction);
	printf("free(ind_cell_multipl);\n");
	free(ind_cell_multipl);
	printf("free(mass1);\n");
	free(mass1);
	printf("free(ind_multipl);\n");
	free(ind_multipl);
	printf("free(snow_region);\n");
	free(snow_region);
	printf("free(bl_cond);\n");
	free(bl_cond);
	printf("free(ind);\n");
	free(ind);
	printf("free(mass);\n");
	free(mass);
	printf("free(A_stencil);\n");
	free(A_stencil);
	printf("free(Aelem_csc);\n");
	free(Aelem_csc);
	printf("free(Ajptr_csc);\n");
	free(Ajptr_csc);
	printf("free(Aiptr_csc);\n");
	free(Aiptr_csc);
	return 0;
}

double density(int i, int j, int k)
{
	double x;
	x = (phase_fraction[PHASE_IND(i, j, k)] * density_snow + (1 - phase_fraction[PHASE_IND(i, j, k)]) * density_air);
	return x;
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

/*
numbering faces of cells is so
0 face - from end of x-axis
1 face - from beginning of x-axis
2 face - from end of y-axis
3 face - from beginning of y-axis
4 face - from end of z-axis
5 face - from beginning of z-axis
*/

int set_arrays(void)
{
	int i, j, k, l, m;
	phase_fraction_len = n_cells_multipl * nz;
	if ((phase_fraction = (double *) malloc(n_cells_multipl * nz * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) phase_fraction, 0, n_cells_multipl * nz * sizeof(double));
	printf("1\n");
	velocity_len = 3 * n_cells_multipl * nz;
	if ((velocity = (double *) malloc(3 * n_cells_multipl * nz * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) velocity, 0, 3 * n_cells_multipl * nz * sizeof(double));
	printf("2\n");
	pressure_len = n_cells_multipl * nz;
	if ((pressure = (double *) malloc(n_cells_multipl * nz * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) pressure, 0, n_cells_multipl * nz * sizeof(double));
	printf("3\n");
	//printf("Setting arrays\n");
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nrows - 1; i++) {
			for (l = 0; l < kx; l++) {
				for (j = 0; j < ncols - 1; j++) {
					for (m = 0; m < ky; m++) {
						if (bl_cond[i * (ncols - 1) + j] != -1) {
							//printf("k = %d, i = %d, l = %d, j = %d, m = %d\n", k, i, l, j, m);
							//printf("snow_region[0] = %d", snow_region[0]);
							//printf("snow_region[i * ncols + j] = %d", snow_region[i * ncols + j]);
							if (snow_region[i * ncols + j] == 1) {
								if ((double) (k + 1) * cellsize >= depth) {
									phase_fraction[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] = 1;
									pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										pressure_atmosphere + density_snow * g[2] * (depth - (double) (k * cellsize));
								} else if (((double) (k + 1) * cellsize < depth) && ((double) k * cellsize) > depth) {
									phase_fraction[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										(depth - (double) k * cellsize) / (double) cellsize;
									pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										pressure_atmosphere + density_snow * g[2] * (depth - (double) (k * cellsize));
								} else {
									phase_fraction[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] = 0;
									pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										pressure_atmosphere;
								}
							} else {
								pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
									pressure_atmosphere;
							}
						}
					}
				}
			}
		}
	}
	printf("Set phase_fraction array\n");
	volume_len = n_cells_multipl;
	if ((volume = (double *) malloc(n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) volume, 0, n_cells_multipl * sizeof(double));
	printf("4\n");
	normal_len = 6 * 3 * n_cells_multipl;
	if ((normal = (double *) malloc(6 * 3 * n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	printf("4.5\n");
	memset((void *) normal, 0, 6 * 3 * n_cells_multipl * sizeof(double));
	printf("5\n");
	area_len = 6 * n_cells_multipl;
	if ((area = (double *) malloc(6 * n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) area, 0, 6 * n_cells_multipl * sizeof(double));
	printf("6\n");
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
				normal[k * 18 + 5 * 3 + 0] = - normal[k * 18 + 0 * 3 + 0];
				normal[k * 18 + 5 * 3 + 1] = - normal[k * 18 + 0 * 3 + 1];
				normal[k * 18 + 5 * 3 + 2] = - normal[k * 18 + 0 * 3 + 2];
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
list of boundary condition modes
	zero_gradient
	no_slip
*/

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

void DDT_density_velocity(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int p;
	for (p = 0; p < 3; p++) {
		A[A_IND(p, i, j, k, p, i, j, k)] += density(i, j, k) / dt;
		A_stencil[A_IND(p, i, j, k, p, i, j, k)] = 1;
		B[B_IND(p, i, j, k)] += density(i, j, k) * velocity[VEL_IND(p, i, j, k)] / dt;
	}
	return;
}

double velocity_on_face(int p, int i, int j, int k, int s)
{
	//printf("velocity_on_face\n");
	if (check_for_corrupt_cell(i, j, k)) return 0;
	int a, b;
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i + 1, j, k)]) / 2;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i - 1, j, k)]) / 2;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i, j + 1, k)]) / 2;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i, j - 1, k)]) / 2;
		case 4:
			if (k + 1 >= nz)
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i, j, k + 1)]) / 2;
		case 5:
			if (k - 1 < 0)
				return 0;
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i, j, k - 1)]) / 2;
	}
}

double density_on_face(int i, int j, int k, int s)
{
	//printf("density_on_face\n");
	if (check_for_corrupt_cell(i, j, k)) return 0;
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i + 1, j, k)) / 2;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i - 1, j, k)) / 2;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i, j + 1, k)) / 2;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i, j - 1, k)) / 2;
		case 4:
			if (k + 1 >= nz)
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i, j, k + 1)) / 2;
		case 5:
			if (k - 1 < 0)
				return 0;
			else
				return (density(i, j, k) + density(i, j, k - 1)) / 2;
	}
}

void DIV_density_velocity_velocity_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s, p, pr;
	for (pr = 0; pr < 3; pr++) {
		for (s = 0; s < 6; s++) {
			for (p = 0; p < 3; p++) {
				B[B_IND(pr, i, j, k)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
					velocity_on_face(p, i, j, k, s) * normal[NORMAL_IND(p, i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)]);
			}
		}
	}
	return;
}

int A_IND_S(int p_eq, int i_eq, int j_eq, int k_eq, int p_el, int i_el, int j_el, int k_el, int s)
{
	switch (s) {
		case 0:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el + 1, j_el, k_el);
		case 1:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el - 1, j_el, k_el);
		case 2:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el, j_el + 1, k_el);
		case 3:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el, j_el - 1, k_el);
		case 4:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el, j_el, k_el + 1);
		case 5:
			return A_IND(p_eq, i_eq, j_eq, k_eq, p_el, i_el, j_el, k_el - 1);
	}
}

int A_IND_S_SWITCH(int i, int j, int k, int s)
{
	if (check_for_corrupt_cell(i, j, k)) return -1;
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
}

void DIV_density_velocity_velocity_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s, p, pr;
	for (pr = 0; pr < 3; pr++) {
		for (s = 0; s < 6; s++) {
			for (p = 0; p < 3; p++) {
				A[A_IND(pr, i, j, k, p, i, j, k)] += area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
					normal[NORMAL_IND(p, i, j, k, s)] / ( 2 * 2 * volume[VOLUME_IND(i, j, k)]);
				A_stencil[A_IND(pr, i, j, k, p, i, j, k)] = 1;
				if (A_IND_S_SWITCH(i, j, k, s) == 1)
				{
					A[A_IND_S(pr, i, j, k, p, i, j, k, s)] += area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
						normal[NORMAL_IND(p, i, j, k, s)] / (2 * 2 * volume[VOLUME_IND(i, j, k)]);
					A_stencil[A_IND_S(pr, i, j, k, p, i, j, k, s)] = 1;
				}
				if (A_IND_S_SWITCH(i, j, k, s) == 0)
				{
					A[A_IND(pr, i, j, k, p, i, j, k)] += area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
						normal[NORMAL_IND(p, i, j, k, s)] / (2 * 2 * volume[VOLUME_IND(i, j, k)]);
					A_stencil[A_IND(pr, i, j, k, p, i, j, k)] = 1;
				}
			}
		}
	}
	return;
}

void DIV_density_velocity_velocity_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	DIV(i, j, k, density_velocity_velocity, half_forward_euler);
	DIV(i, j, k, density_velocity_velocity, half_backward_euler);
	return;
}

void GRAD_pressure_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	GRAD(i, j, k, pressure, half_forward_euler);
	GRAD(i, j, k, pressure, half_backward_euler);
	return;
}

void GRAD_pressure_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1)) {
		B[B_IND(0, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure_atmosphere) / (2 * dx);
	} else {
		B[B_IND(0, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure[PRESS_IND(i - 1, j, k)]) / (2 * dx);
	}
	if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1)) {
		B[B_IND(1, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure_atmosphere) / (2 * dy);
	} else {
		B[B_IND(1, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure[PRESS_IND(i, j - 1, k)]) / (2 * dy);
	}
	if (k - 1 >= 0) {
		B[B_IND(2, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure[PRESS_IND(i, j, k - 1)]) / (2 * dz);
	}
	return;
}

void GRAD_pressure_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1)) {
		A[A_IND(0, i, j, k, 4, i, j, k)] += 1 / (2 * dx);
		A_stencil[A_IND(0, i, j, k, 4, i, j, k)] = 1;
		B[B_IND(0, i, j, k)] += pressure_atmosphere / (2 * dx);
	} else {
		A[A_IND(0, i, j, k, 4, i, j, k)] += 1 / (2 * dx);
		A_stencil[A_IND(0, i, j, k, 4, i, j, k)] = 1;
		A[A_IND(0, i, j, k, 4, i - 1, j, k)] -= 1 / (2 * dx);
		A_stencil[A_IND(0, i, j, k, 4, i - 1, j, k)] = 1;
	}
	if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1)) {
		A[A_IND(1, i, j, k, 4, i, j, k)] += 1 / (2 * dy);
		A_stencil[A_IND(1, i, j, k, 4, i, j, k)] = 1;
		B[B_IND(1, i, j, k)] += pressure_atmosphere / (2 * dy);
	} else {
		A[A_IND(1, i, j, k, 4, i, j, k)] += 1 / (2 * dy);
		A_stencil[A_IND(1, i, j, k, 4, i, j, k)] = 1;
		A[A_IND(1, i, j, k, 4, i, j - 1, k)] -= 1 / (2 * dy);
		A_stencil[A_IND(1, i, j, k, 4, i, j - 1, k)] = 1;
	}
	if (k - 1 >= 0) {
		A[A_IND(2, i, j, k, 4, i, j, k)] += 1 / (2 * dz);
		A_stencil[A_IND(2, i, j, k, 4, i, j, k)] = 1;
		A[A_IND(2, i, j, k, 4, i, j, k - 1)] -= 1 / (2 * dz);
		A_stencil[A_IND(2, i, j, k, 4, i, j, k - 1)] = 1;
	}
	return;
}

void VECT_gravity_force_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int p;
	for (p = 0; p < 3; p++) {
		B[B_IND(p, i, j, k)] += density_snow * g[p] * phase_fraction[PHASE_IND(i, j, k)] / 2 + density_air * g[p] * (1 - phase_fraction[PHASE_IND(i, j, k)]) / 2;
	}
	return;
}

void VECT_gravity_force_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int p;
	for (p = 0; p < 3; p++) {
		A[A_IND(p, i, j, k, 3, i, j, k)] -= density_snow * g[p] / 2 - density_air * g[p] / 2;
		A_stencil[A_IND(p, i, j, k, 3, i, j, k)] = 1;
		B[B_IND(p, i, j, k)] += density_air * g[p] / 2;
	}
	return;
}

void VECT_gravity_force_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	VECT(i, j, k, gravity_force, half_forward_euler);
	VECT(i, j, k, gravity_force, half_backward_euler);
	return;
}

double strain_rate_on_face(int i, int j, int k, int s, int m, int n)
{
	if (check_for_corrupt_cell(i, j, k)) return 0;
	switch (m + n + m * n) {
		case 0:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return 0;
			else
				return (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		case 1:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
					return 0;
				else
					return 0.5 * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / dx;
			else
				if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
					return 0.5 * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j - 1, k, s)) / dy;
				else
					return 0.5 * ((velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j - 1, k, s)) / dy + (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / dx);
		case 2:
			if (k - 1 < 0)
				if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
					return 0.5 * velocity_on_face(0, i, j, k, s) / dz;
				else
					return 0.5 * (velocity_on_face(0, i, j, k, s) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / dx);
			else
				if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
					return 0.5 * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j, k - 1, s)) / dz;
				else
					return 0.5 * ((velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j, k - 1, s)) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / dx);
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return 0;
			else
				return (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j - 1, k, s)) / dx;
		case 5:
			if (k - 1 < 0)
				if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
					return 0.5 * velocity_on_face(1, i, j, k, s) / dz;
				else
					return 0.5 * (velocity_on_face(1, i, j, k, s) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j - 1, k, s)) / dy);
			else
				if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
					return 0.5 * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j, k - 1, s)) / dz;
				else
					return 0.5 * ((velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j, k - 1, s)) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j - 1, k, s)) / dy);
		case 8:
			if (k <= 0)
				return velocity_on_face(2, i, j, k, s) / dx;
			else
				return (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j, k - 1, s)) / dx;
	}
}

double shear_rate_on_face(int i, int j, int k, int s)
{
	if (check_for_corrupt_cell(i, j, k)) return 0;
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
	if (check_for_corrupt_cell(i, j, k)) return 1;
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i + 1, j, k)]) / 2;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i - 1, j, k)]) / 2;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i, j + 1, k)]) / 2;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i, j - 1, k)]) / 2;
		case 4:
			if (k + 1 >= nz)
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i, j, k + 1)]) / 2;
		case 5:
			if (k - 1 < 0)
				return 1;
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i, j, k - 1)]) / 2;
	}
}

double effective_viscosity_on_face(int i, int j, int k, int s)
{
	if (check_for_corrupt_cell(i, j, k)) return 0;
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

void DIV_shear_stress_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s, p, m;
	for (p = 0; p < 3; p++) {
		for (s = 0; s < 6; s++) {
			for (m = 0; m < 3; m++) {
				B[B_IND(p, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * 2 * effective_viscosity_on_face(i, j, k, s) * strain_rate_on_face(i, j, k, s, p, m) * normal[NORMAL_IND(m, i, j, k, s)];
			}
		}
	}
	return;
}

double pressure_on_face(int i, int j, int k, int s)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i + 1, j, k)]) / 2;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i - 1, j, k)]) / 2;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i, j + 1, k)]) / 2;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i, j - 1, k)]) / 2;
		case 4:
			if (k + 1 >= nz)
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i, j, k + 1)]) / 2;
		case 5:
			if (k - 1 < 0)
				return pressure[PRESS_IND(i, j, k)];
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i, j, k - 1)]) / 2;
	}
}

void DIV_grad_pressure_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	DIV(i, j, k, grad_pressure, half_forward_euler);
	DIV(i, j, k, grad_pressure, half_backward_euler);
	return;
}

void DIV_grad_pressure_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_atmosphere) * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		else
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i - 1, j, k, s)) * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_atmosphere) * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		else
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i, j - 1, k, s)) * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		if (k > 0)
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i, j, k - 1, s)) * normal[NORMAL_IND(2, i, j, k, s)] / dz;
	}
	return;
}

void DIV_grad_pressure_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1)) {
			A[A_IND(4, i, j, k, 4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k)] = 1;
			if (A_IND_S_SWITCH(i, j, k, s)) {
				A[A_IND_S(4, i, j, k, 4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		} else {
			A[A_IND(4, i, j, k, 4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k)] = 1;
			A[A_IND(4, i, j, k, 4, i - 1, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			A_stencil[A_IND(4, i, j, k, 4, i - 1, j, k)] = 1;
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			if (A_IND_S_SWITCH(i - 1, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i - 1, j, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
				A_stencil[A_IND_S(4, i, j, k, 4, i - 1, j, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		}
		if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1)) {
			A[A_IND(4, i, j, k, 4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k)] = 1;
			if (A_IND_S_SWITCH(i, j, k, s)) {
				A[A_IND_S(4, i, j, k, 4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		} else {
			A[A_IND(4, i, j, k, 4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k)] = 1;
			A[A_IND(4, i, j, k, 4, i, j - 1, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			A_stencil[A_IND(4, i, j, k, 4, i, j - 1, k)] = 1;
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			if (A_IND_S_SWITCH(i, j - 1, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i, j - 1, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j - 1, k, s)] = 1;
			} else
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		}
		if (k > 0) {
			A[A_IND(4, i, j, k, 4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k)] = 1;
			A[A_IND(4, i, j, k, 4, i, j, k - 1)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			A_stencil[A_IND(4, i, j, k, 4, i, j, k - 1)] = 1;
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k, s)] = 1;
			} else {
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			}
			if (A_IND_S_SWITCH(i, j, k - 1, s) == 1) {
				A[A_IND_S(4, i, j, k, 4, i, j, k - 1, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				A_stencil[A_IND_S(4, i, j, k, 4, i, j, k - 1, s)] = 1;
			} else {
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			}
		}
	}
	return;
}

void DIV_div_density_velocity_velocity_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	DIV(i, j, k, div_density_velocity_velocity, half_forward_euler);
	DIV(i, j, k, div_density_velocity_velocity, half_backward_euler);
	return;
}

void DIV_div_density_velocity_velocity_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		if (!((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1)))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		if (!((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1)))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j - 1, k, s)) / dy;
		if (k <= 0)
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * velocity_on_face(2, i, j, k, s) / dz;
		else	
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j, k - 1, s)) / dz;
	}
	return;
}

void DIV_div_density_velocity_velocity_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		if (!((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))) {
			A[A_IND(4, i, j, k, 0, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			A_stencil[A_IND(4, i, j, k, 0, i, j, k)] = 1; 
			A[A_IND(4, i, j, k, 0, i - 1, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			A_stencil[A_IND(4, i, j, k, 0, i - 1, j, k)] = 1; 
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 0, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
				A_stencil[A_IND_S(4, i, j, k, 0, i, j, k, s)] = 1;
			} else {
				A[A_IND(4, i, j, k, 0, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
				A_stencil[A_IND(4, i, j, k, 0, i, j, k)] = 1;
			}
			if (A_IND_S_SWITCH(i - 1, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 0, i - 1, j, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
				A_stencil[A_IND_S(4, i, j, k, 0, i - 1, j, k, s)] = 1;
			} else {
				A[A_IND(4, i, j, k, 0, i - 1, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
				A_stencil[A_IND(4, i, j, k, 0, i - 1, j, k)] = 1;
			}
		}
		if (!((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))) {
			A[A_IND(4, i, j, k, 1, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			A_stencil[A_IND(4, i, j, k, 1, i, j, k)] = 1; 
			A[A_IND(4, i, j, k, 1, i, j - 1, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			A_stencil[A_IND(4, i, j, k, 1, i, j - 1, k)] = 1; 
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 1, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
				A_stencil[A_IND_S(4, i, j, k, 1, i, j, k, s)] = 1; 
			} else {
				A[A_IND(4, i, j, k, 1, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
				A_stencil[A_IND(4, i, j, k, 1, i, j, k)] = 1;
			}
			if (A_IND_S_SWITCH(i, j - 1, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 1, i, j - 1, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
				A_stencil[A_IND_S(4, i, j, k, 1, i, j - 1, k, s)] = 1; 
			} else {
				A[A_IND(4, i, j, k, 1, i, j - 1, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
				A_stencil[A_IND(4, i, j, k, 1, i, j - 1, k)] = 1;
			}
		}
		if (k == 0) {
			A[A_IND(4, i, j, k, 2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			A_stencil[A_IND(4, i, j, k, 2, i, j, k)] = 1; 
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 2, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND_S(4, i, j, k, 2, i, j, k, s)] = 1;
			} else {
				A[A_IND(4, i, j, k, 2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND(4, i, j, k, 2, i, j, k)] = 1;
			}
		}
		if (k > 0) {
			A[A_IND(4, i, j, k, 2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			A_stencil[A_IND(4, i, j, k, 2, i, j, k)] = 1; 
			A[A_IND(4, i, j, k, 2, i, j, k - 1)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			A_stencil[A_IND(4, i, j, k, 2, i, j, k - 1)] = 1; 
			if (A_IND_S_SWITCH(i, j, k, s) == 1) {
				A[A_IND_S(4, i, j, k, 2, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND_S(4, i, j, k, 2, i, j, k, s)] = 1; 
			} else {
				A[A_IND(4, i, j, k, 2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND(4, i, j, k, 2, i, j, k)] = 1;
			}
			if (A_IND_S_SWITCH(i, j, k - 1, s) == 1) {
				A[A_IND_S(4, i, j, k, 2, i, j, k - 1, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND_S(4, i, j, k, 2, i, j, k - 1, s)] = 1; 
			} else {
				A[A_IND(4, i, j, k, 2, i, j, k - 1)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
				A_stencil[A_IND(4, i, j, k, 2, i, j, k - 1)] = 1;
			}
		}
	}
	return;
}

void DDT_density_snow_volume_fraction(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	A[A_IND(3, i, j, k, 3, i, j, k)] += density(i, j, k) / dt;
	A_stencil[A_IND(3, i, j, k, 3, i, j, k)] = 1;
	B[B_IND(3, i, j, k)] += density(i, j, k) * phase_fraction[PHASE_IND(i, j, k)] / dt;
	return;
}

void DIV_density_snow_volume_fraction_velocity_crank_nikolson(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	DIV(i, j, k, density_snow_volume_fraction_velocity, half_forward_euler);
	DIV(i, j, k, density_snow_volume_fraction_velocity, half_backward_euler);
	return;
}

void DIV_density_snow_volume_fraction_velocity_half_forward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		B[B_IND(3, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * phase_fraction_on_face(i, j, k, s);
	}
	return;
}

void DIV_density_snow_volume_fraction_velocity_half_backward_euler(int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return;
	int s;
	for (s = 0; s < 6; s++) {
		A[A_IND(3, i, j, k, 3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 4);
		A_stencil[A_IND(3, i, j, k, 3, i, j, k)] = 1;
		if (A_IND_S_SWITCH(i, j, k, s) == 1) {
			A[A_IND_S(3, i, j, k, 3, i, j, k, s)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 4);
			A_stencil[A_IND_S(3, i, j, k, 3, i, j, k, s)] = 1;
		}
		if (A_IND_S_SWITCH(i, j, k, s) == 0) {
			A[A_IND(3, i, j, k, 3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 4);
			A_stencil[A_IND(3, i, j, k, 3, i, j, k)] = 1;
		}
	}
	return;
}

int create_Ab(void)
{
	num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	system_dimension = n_cells_multipl * nz * num_parameters;
	A_len = system_dimension * system_dimension;
	if ((A = (double *) malloc(system_dimension * system_dimension * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) A, 0, system_dimension * system_dimension * sizeof(double));
	if ((A_stencil = (int *) malloc(system_dimension * system_dimension * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) A_stencil, 0, system_dimension * system_dimension * sizeof(int));
	B_len = system_dimension;
	if ((B = (double *) malloc(system_dimension * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) B, 0, system_dimension * sizeof(double));
	int i, j, k;
	set_arrays();
	if (print_vtk(0) == 1)
		printf("Error printing vtk file\n");
	printf("we came here?\n");
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1) {
					DDT(i, j, k, density_velocity);
					DIV(i, j, k, density_velocity_velocity, crank_nikolson);
					VECT(i, j, k, gravity_force, crank_nikolson);
					GRAD(i, j, k, pressure, crank_nikolson);
					DIV(i, j, k, shear_stress, forward_euler);
					DIV(i, j, k, grad_pressure, crank_nikolson);
					DIV(i, j, k, div_density_velocity_velocity, crank_nikolson);
					DDT(i, j, k, density_snow_volume_fraction);
					DIV(i, j, k, density_snow_volume_fraction_velocity, crank_nikolson);
				}
			}
		}
	}
	printf_Ab();
	solve_matrix();
	if (print_vtk(1) == 1)
		printf("Error printing vtk file\n");
	if (A_csr_func() == 1) return 1;
	return 0;
}

void printf_Ab(void)
{
	FILE *f = fopen("A.txt","w");
	int i, j;
	for (i = 0; i < system_dimension; i++) {
		for (j = 0; j < system_dimension; j++) {
			fprintf(f, "%010lf\t", A[i * system_dimension + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen("b.txt","w");
	for (i = 0; i < system_dimension; i++) {
		fprintf(f, "%010lf\n", B[i]);
	}
	fclose(f);
	f = fopen("A_portrait.txt", "w");
	for (i = 0; i < system_dimension; i++) {
		for (j = 0; j < system_dimension; j++) {
			fprintf(f, "%d\t%d\t", i, j);
			if (A[i * system_dimension + j] == 0)
				fprintf(f, "0\n");
			else
				fprintf(f, "1\n");
		}
	}
	fclose(f);
	return;
}

int A_csc_func(void)
{
	int i, j, k = 0, l;
	non_zero_elem = 0;
	for (i = 0; i < system_dimension * system_dimension; i++)
		if (A_stencil[i] != 0) non_zero_elem++;
	if ((Aelem_csc = (double *) malloc(non_zero_elem * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Aelem_csc, 0, non_zero_elem * sizeof(double));
	if ((Aiptr_csc = (int *) malloc(non_zero_elem * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Aiptr_csc, 0, non_zero_elem * sizeof(double));
	if ((Ajptr_csc = (int *) malloc((system_dimension + 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Ajptr_csc, 0, (system_dimension + 1) * sizeof(double));
	for (j = 0; j < system_dimension; j++) {
		l = 0;
		for (i = 0; i < system_dimension; i++) {
			if (A_stencil[i * system_dimension + j] != 0) {
				Aelem_csc[k] = A[i * system_dimension + j];
				Aiptr_csc[k] = i;
				if (l == 0) {
					Ajptr_csc[j] = k;
					l = 1;
				}
				k++;
			}
		}
		if (l == 0)
			Ajptr_csc[j] = Ajptr_csc[j - 1];
	}
	Ajptr_csc[j] = k;
//	print_A_csr();
	return 0;
}

int A_csr_func(void)
{
	int i, j, k = 0, l;
	non_zero_elem = 0;
	for (i = 0; i < system_dimension * system_dimension; i++)
		if (A_stencil[i] != 0) non_zero_elem++;
	if ((Aelem_csr = (double *) malloc(non_zero_elem * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Aelem_csr, 0, non_zero_elem * sizeof(double));
	if ((Ajptr_csr = (int *) malloc(non_zero_elem * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Ajptr_csr, 0, non_zero_elem * sizeof(double));
	if ((Aiptr_csr = (int *) malloc((system_dimension + 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset(Aiptr_csr, 0, (system_dimension + 1) * sizeof(double));
	for (i = 0; i < system_dimension; i++) {
		l = 0;
		for (j = 0; j < system_dimension; j++) {
			if (A_stencil[i * system_dimension + j] != 0) {
				Aelem_csr[k] = A[i * system_dimension + j];
				Ajptr_csr[k] = j;
				if (l == 0) {
					Aiptr_csr[i] = k;
					l = 1;
				}
				k++;
			}
		}
		if (l == 0)
			Aiptr_csr[i] = Aiptr_csr[i - 1];
	}
	Aiptr_csr[i] = k;
	print_A_csr();
	return 0;
}

void print_A_csr(void)
{
	printf("printf_A_csr\n");
	FILE *f;
	int i;
	printf("file declared\n");
	if ((f = fopen("A_csr_elem.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	printf("file opened\n");
	printf("8\n");
	for (i = 0; i < non_zero_elem; i++) {
		fprintf(f, "%010lf\n", Aelem_csr[i]);
	}
	fclose(f);
	printf("9\n");
	f = fopen("A_csr_jptr.txt","w");
	for (i = 0; i < non_zero_elem; i++) {
		fprintf(f, "%010d\n", Ajptr_csr[i]);
	}
	fclose(f);
	printf("10\n");
	f = fopen("A_csr_iptr.txt","w");
	for (i = 0; i < system_dimension; i++) {
		fprintf(f, "%010d\n", Aiptr_csr[i]);
	}
	fclose(f);
	return;
}

int solve_matrix(void)
{
	SuperMatrix A_csc;
	NCformat *Astore;
//	double   *a;
//	int      *asub, *xa;
	int      *perm_c; /* column permutation vector */
	int      *perm_r; /* row permutations from partial pivoting */
	SuperMatrix L;      /* factor L */
	SCformat *Lstore;
	SuperMatrix U;      /* factor U */
	NCformat *Ustore;
	SuperMatrix B_csc;
//	int      nrhs, ldx, info, m, n, nnz;
	int      info, i, nrhs;
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
	dCreate_CompCol_Matrix(&A_csc, system_dimension, system_dimension, non_zero_elem, Aelem_csc, Aiptr_csc, Ajptr_csc, SLU_NC, SLU_D, SLU_GE);
	Astore = A_csc.Store;
	printf("Dimension %dx%d; # nonzeros %d\n", A_csc.nrow, A_csc.ncol, Astore->nnz);
	nrhs   = 1;
//	if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B_csc, system_dimension, nrhs, B, system_dimension, SLU_DN, SLU_D, SLU_GE);
//	xact = doubleMalloc(n * nrhs);
//	ldx = n;
//	dGenXtrue(n, nrhs, xact, ldx);
//	dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
	if ( !(perm_c = intMalloc(system_dimension)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(perm_r = intMalloc(system_dimension)) ) ABORT("Malloc fails for perm_r[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, &A_csc, perm_c, perm_r, &L, &U, &B_csc, &stat, &info);
	if ( info == 0 ) {
		/* This is how you could access the solution matrix. */
		double *sol = (double*) ((DNformat*) B_csc.Store)->nzval; 
		/* Compute the infinity norm of the error. */
//		dinf_norm_error(nrhs, &B, xact);
		dinf_norm_error(nrhs, &B_csc, sol);
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
	/*
	for (i = 0; i < n_cells_multipl * nz; i++) {
		velocity[i * 3 + 0] = sol[i * num_parameters + 0];
		velocity[i * 3 + 1] = sol[i * num_parameters + 1];
		velocity[i * 3 + 3] = sol[i * num_parameters + 2];
		phase_fraction[i] = sol[i * num_parameters + 3];
		pressure[i] = sol[i * num_parameters + 4];
	}
	*/
	
	StatFree(&stat);
//	SUPERLU_FREE (rhs);
//	SUPERLU_FREE (xact);
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);
	Destroy_CompCol_Matrix(&A_csc);
	Destroy_SuperMatrix_Store(&B_csc);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Exit main()");
#endif
	return 0;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S must be declared\n\n");
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
	printf("-h\tdisplay usage\n");
}

int main(int argc, char **argv)
{
	int opt = 0;
	g[0] = g[1] = 0;
	g[2] = 9,81;
	static const char *optString = "m:r:H:D:x:y:z:a:s:p:v:k:i:l:S:h?";
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
			case 'h':
			case '?':
				display_usage();
				return 1;
		}
	}
	if (argc != 31) {
		printf("Not enouth arguments\n");
		goto error;
	}
	printf("Parameters entered correctly.\n");

	if (read_asc_and_declare_variables() == 1) goto error;
	printf("Files of map and region have been processed\n");
	if (do_interpolation() == 1) goto error;
	printf("The interpolation have been done\n");
	dt = 0.01;//we need to set dt!!!
	printf("Creating Ax = b equation\n");
	if (create_Ab() == 1) goto error;
	printf("Ax = b equation well done\n");
	printf("Free massives\n");
	if (free_massives() == 1) goto error;

	printf("Calculations finished successfully\n");
	return 0;
error:
	printf("Error\n");
	display_usage();
	return 1;
}
