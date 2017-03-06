#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

char *map_name;
char *region_map_name;
float hight;
int kx, ky, kz;
int nx, ny, nz;
int n_bl_x, n_bl_y, n_bl_z;
int ncols;
int nrows;
float cellsize;
float *mass;
int *ind;
int *bl_cond;
int n_points, n_cells;
int n_points_multipl, n_cells_multipl;
int *ind_multipl, *ind_cell_multipl;
float interpolation, interpolation_poli;
float *mass1;
float nodata_value;
int *snow_region;
float depth;
float density_snow, density_air;
float *phase_fraction;
float *velocity;
float *pressure;
float pressure_atmosphere;
float g[3];
float k_viscosity_snow, k_viscosity_air;
float *normal, *volume, *area;
float dt, dx, dy, dz;
float *A, *B;
int num_parameters;
float shear_rate_0,limiting_viscosity_snow, flow_index, yield_stress;

int read_asc_and_declare_variables(void);
void annihilate_array(void *a, int size_bites);
int do_interpolation(void);
int print_vtk(void);
int free_massives(void);
float density(int i, int j, int k);
void print_mass(void);
int set_arrays(void);
int conformity(int i, int j);
void DDT_density_velocity(int i, int j, int k);
float velocity_on_face(int p, int i, int j, int k, int s);
float density_on_face(int i, int j, int k, int s);
void DIV_density_velocity_velocity_half_forward_euler(int i, int j, int k);
int A_IND_S(int p, int i, int j, int k, int s);
int A_IND_S_SWITCH(int i, int j, int k, int s);
void DIV_density_velocity_velocity_half_backward_euler(int i, int j, int k);
void DIV_density_velocity_velocity_crank_nikolson(int i, int j, int k);
void GRAD_pressure_crank_nikolson(int i, int j, int k);
void GRAD_pressure_half_forward_euler(int i, int j, int k);
void GRAD_pressure_half_backward_euler(int i, int j, int k);
void VECT_gravity_force_half_forward_euler(int i, int j, int k);
void VECT_gravity_force_half_backward_euler(int i, int j, int k);
void VECT_gravity_force_crank_nikolson(int i, int j, int k);
float strain_rate_on_face(int i, int j, int k, int s, int m, int n);
float shear_rate_on_face(int i, int j, int k, int s);
float phase_fraction_on_face(int i, int j, int k, int s);
float effective_viscosity_on_face(int i, int j, int k, int s);
void DIV_shear_stress_forward_euler(int i, int j, int k);
float pressure_on_face(int i, int j, int k, int s);
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
int A_IND(int p, int i, int j, int k);
int B_IND(int p, int i, int j, int k);
int PRESS_IND(int i, int j, int k);
int VEL_IND(int p, int i, int j, int k);
int AREA_IND(int i, int j, int k, int s);
int NORMAL_IND(int p, int i, int j, int k, int s);
int VOLUME_IND(int i, int j, int k);
int PHASE_IND(int i, int j, int k);

/*
#define A_IND(p, i, j, k) \
	((num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p) * \
	n_cells_multipl * nz * num_parameters + \
	num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p)

#define B_IND(p, i, j, k) \
	(num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p)

#define PRESS_IND(i, j, k) \
	(k * n_cells_multipl + ind_cell_multipl[i * ny + j])

#define VEL_IND(p, i , j, k) \
	(3 * PRESS_IND(i, j, k) + p)

#define AREA_IND(i, j, k, s) \
	(6 * ind_cell_multipl[i * ny + j] + s)

#define NORMAL_IND(p, i, j, k, s) \
	(18 * ind_cell_multipl[i * ny + j] + s * 3 + p)

#define VOLUME_IND(i, j, k) \
	(ind_cell_multipl[i * ny + j])

#define PHASE_IND(i, j, k) \
	(k * n_cells_multipl + ind_cell_multipl[i * ny + j])
*/

#define DDT(i, j, k, object) DDT_##object(i, j, k);
#define DIV(i, j, k, object, numerical_scheme) DIV_##object##_##numerical_scheme(i, j, k);
#define BOUNDARY_CONDITION(p, i, j, k, s, object, mode, numerical_scheme) BOUNDARY_CONDITION_##object##_##mode##_##numerical_scheme(p, i, j, k, s)
#define DDX(p, i, j, k, s, mode, numerical_scheme) DDX_##mode##_##numerical_scheme(p, i, j, k, s)
#define GRAD(i, j, k, mode, numerical_scheme) GRAD_##mode##_##numerical_scheme(i, j, k)
#define VECT(i, j, k, mode, numerical_scheme) VECT_##mode##_##numerical_scheme(i, j, k)
#define DDY(p, i, j, k, s, mode, numerical_scheme) DDY_##mode##_##numerical_scheme(p, i, j, k, s)
#define DDZ(p, i, j, k, s, mode, numerical_scheme) DDZ_##mode##_##numerical_scheme(p, i, j, k, s)

int A_IND(int p, int i, int j, int k)
{
	//printf("A_IND p = %d, i = %d, j = %d, k = %d\n", p, i, j, k);
	return ((num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p) *
		n_cells_multipl * nz * num_parameters + num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p);
}

int B_IND(int p, int i, int j, int k)
{
	return (num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p);
}

int PRESS_IND(int i, int j, int k)
{
	return (k * n_cells_multipl + ind_cell_multipl[i * ny + j]);
}

int VEL_IND(int p, int i, int j, int k)
{
	return (3 * PRESS_IND(i, j, k) + p);
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

int PHASE_IND(int i, int j, int k)
{
	return (k * n_cells_multipl + ind_cell_multipl[i * ny + j]);
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

	float xllcorner;
	if ((err = fscanf(f, "%s %f", str, &xllcorner)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "xllcorner") != 0) {
		printf("Error file of map: no \"xllcorner\" tag\n");
		return 1;
	}
	
	float yllcorner;
	if ((err = fscanf(f, "%s %f", str, &yllcorner)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "yllcorner") != 0) {
		printf("Error file of map: no \"yllcorner\" tag\n");
		return 1;
	}
	
	if ((err = fscanf(f, "%s %f", str, &cellsize)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "cellsize") != 0) {
		printf("Error file of map: no \"cellsize\" tag\n");
		return 1;
	}

	if ((err = fscanf(f, "%s %f", str, &nodata_value)) == EOF) {
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

	float xllcorner1;
	if ((err = fscanf(f1, "%s %f", str, &xllcorner1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "xllcorner") != 0) {
		printf("Error file of region: no \"xllcorner\" tag\n");
		return 1;
	}
	
	float yllcorner1;
	if ((err = fscanf(f1, "%s %f", str, &yllcorner1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "yllcorner") != 0) {
		printf("Error file of region: no \"yllcorner\" tag\n");
		return 1;
	}
	
	float cellsize1;
	if ((err = fscanf(f1, "%s %f", str, &cellsize1)) == EOF) {
		printf("Error file of region\n");
		return 1;
	}
	if (strcmp(str, "cellsize") != 0) {
		printf("Error file of region: no \"cellsize\" tag\n");
		return 1;
	}

	float nodata_value1;
	if ((err = fscanf(f1, "%s %f", str, &nodata_value1)) == EOF) {
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

	if ((mass = (float *) malloc(ncols * nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((ind = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((bl_cond = (int *) malloc((ncols - 1) * (nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	for (i = 0; i < ncols * nrows; i++) {
		if (err = fscanf(f, "%f", &mass[i]) == EOF) {
			printf("Error file of map\n");
			return 1;
		}
		ind[i] = -1;
		if (i < (ncols - 1) * (nrows - 1))
			bl_cond[i] = -1;
	}
	fclose(f);
	remove("map.txt");
	
	float *mass_tmp;
	if ((mass_tmp = (float *) malloc(ncols1 * nrows1 * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	for (i = 0; i < ncols1 * nrows1; i++) {
		if (err = fscanf(f1, "%f", &mass_tmp[i]) == EOF) {
			printf("Error file of region\n");
			return 1;
		}
	}
	remove("regions_map.txt");

	if ((snow_region = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	float xlucorner = xllcorner - cellsize * nrows;
	float ylucorner = yllcorner;
	float xlucorner1 = xllcorner1 - cellsize1 * nrows1;
	float ylucorner1 = yllcorner1;
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

void annihilate_array(void *a, int size_bites)
{
	char *b;
	b = (char *) a;
	int i;
	for (i = 0; i < size_bites; i++)
		b[i] = '0';
	return;
}

int do_interpolation(void)
{
	//printf("Interpolating mesh.\n");
	int i, j, k, l = 0;
	float a, b, d;
	float *c;
	float *e, *f;
	
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

	int count_ind_multipl = 0;
	if ((mass1 = (float *) malloc(((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1) * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
		mass1[i] = nodata_value;
	for (i = 0; i < nrows; i++)
		for (j = 0; j < ncols; j++)
			mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky] = mass[i * ncols + j];

	//printf("mass1 after assignment mass\n");
	//for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
	//	printf("%f\n", mass1[i]);

	if ((c = (float *) malloc(ncols * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (float *) malloc(ncols * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (float *) malloc(ncols * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	int ind_start, ind_finish, flag;
	for (i = 0; i < nrows; i++) {
		flag = 0;
		ind_start = 0;
		ind_finish = ncols;
		for (j = 0; j < ncols; j++) {
			if ((mass[i * ncols + j] != nodata_value) && (flag == 0)) {
				ind_start = j;
				flag = 1;
			}
			if ((mass[i * ncols + j] == nodata_value) && (flag == 1)) {
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
			if ((mass[i * ncols + j - 1] != nodata_value) && (mass[i * ncols + j] != nodata_value)) {
				for (k = 0; k < ky; k++) {
					mass1[i * ((ncols - 1) * ky + 1) * kx + (j - 1) * ky + k] = a + b * ((float) k * cellsize / (float) ky - cellsize) +
						c[j] * pow((float) k * cellsize / (float) ky - cellsize, 2) / 2 +
						d * pow((float) k * cellsize / (float) ky - cellsize, 3) / 6;
				}
			}
		}
		annihilate_array((void *) c, ncols * sizeof(float));
		annihilate_array((void *) e, ncols * sizeof(float));
		annihilate_array((void *) f, ncols * sizeof(float));
	}

	//printf("Interpolation along the y axis have been done\n");
	
	if ((c = (float *) realloc(c, nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (float *) realloc(e, nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (float *) realloc(f, nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	for (j = 0; j < ncols; j++) {
		for (l = 0; l < ky; l++) {
			flag = 0;
			ind_start = 0;
			ind_finish = nrows;
			for (i = 0; i < nrows; i++) {
				if (((mass[i * ncols + j] != nodata_value) && (flag == 0) && (l == 0)) ||
					(((mass[i * ncols + j] != nodata_value) && (mass[i * ncols + j + 1] != nodata_value)) &&
					 	(flag == 0) && (l != 0))) {
							ind_start = i;
							flag = 1;
				}
				if (((mass[i * ncols + j] == nodata_value) && (flag == 1) && (l == 0)) ||
					(((mass[i * ncols + j] == nodata_value) || (mass[i * ncols + j + 1] == nodata_value)) &&
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
				mass1[i * ((ncols - 1) * ky + 1) * kx + j * ky + l] +
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
								a + b * ((float) k * cellsize / (float) kx - cellsize) +
								c[i] * pow((float) k * cellsize / (float) kx - cellsize, 2) / 2 +
								d * pow((float) k * cellsize / (float) kx - cellsize, 3) / 6;
						}
				}
			}
			annihilate_array((void *) c, nrows * sizeof(float));
			annihilate_array((void *) e, nrows * sizeof(float));
			annihilate_array((void *) f, nrows * sizeof(float));
			if (j == ncols - 1) break;
		}
	}

	free(c);
	free(e);
	free(f);

	//printf("Interpolation along the x axes have been done\n");

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

	dx /= (float) kx;
	dy /= (float) ky;
	dz /= (float) kz;
	return 0;
}

int print_vtk(void)
{
	FILE *f = fopen("map.vtk", "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "slope\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

	int i, j, k, a;
	fprintf(f, "POINTS %d float\n", n_points_multipl * ((int) (hight / cellsize) * kz + 1));
	for (k = 0; k <= (int) (hight / cellsize) * kz; k++) {
		for (i = 0; i < (nrows - 1) * kx + 1; i ++) {
			for (j = 0; j < (ncols - 1) * ky + 1; j++) {
				if (ind_multipl[i * ((ncols - 1) * ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * cellsize / (float) kx, j * cellsize / (float) ky,
						mass1[i * ((ncols - 1) * ky + 1) + j] + k * cellsize / (float) kz);
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
	return 0;
}

int free_massives(void)
{
	printf("free(B);\n");
	//free(B);
	printf("free(A);\n");
	free(A);
	printf("free(area);\n");
	free(area);
	printf("free(normal);\n");
	free(normal);
	printf("free(volume);\n");
	free(volume);
	printf("free(pressure);\n");
	//free(pressure);
	printf("free(velocity);\n");
	free(velocity);
	printf("free(phase_fraction);\n");
	//free(phase_fraction);
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
	printf("free(pressure);\n");
	free(pressure);
	printf("free(B);\n");
	free(B);
	return 0;
}

float density(int i, int j, int k)
{
	float x;
	//printf("Calculating density\n");
	//printf("i = %d, j = %d, k = %d\n", i, j, k);
	//printf("phase_fraction[PHASE_IND(i, j, k)] = %f\n", phase_fraction[PHASE_IND(i, j, k)]);
	x = (phase_fraction[PHASE_IND(i, j, k)] * density_snow + (1 - phase_fraction[PHASE_IND(i, j, k)]) * density_air);
	//printf("End\n");
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
	int i, j, k, l, m, n;
	if ((phase_fraction = (float *) malloc(n_cells_multipl * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((velocity = (float *) malloc(3 * n_cells_multipl * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	annihilate_array(velocity, 3 * n_cells_multipl * nz * sizeof(float));
	if ((pressure = (float *) malloc(n_cells_multipl * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	n = 0;
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
								if (((float) (k + 1) * cellsize > depth) && (depth > (float) k * cellsize)) {
									phase_fraction[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										(depth - (float) k * cellsize) / (float) cellsize;
									pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										pressure_atmosphere + density_snow * g[2] * (depth - (float) (k * cellsize));
								}
								else if (depth > (float) (k + 1) * cellsize) {
									phase_fraction[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] = 1;
									pressure[k * n_cells_multipl + ind_cell_multipl[i * (ncols - 1) * ky + j * ky + m]] =
										pressure_atmosphere + density_snow * g[2] * (depth - (float) k * cellsize);
								}
								else {
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
	//printf("Set phase_fraction array\n");
	if ((volume = (float *) malloc(n_cells_multipl * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((normal = (float *) malloc(n_cells_multipl * 6 * 3 * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((area = (float *) malloc(n_cells_multipl * 6 * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	float a[4], b[4], c[4], p, d;
	k = 0;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (ind_cell_multipl[i * ny + j] != -1) {
				a[1] = cellsize / (float) kx;
				a[2] = 0;
				a[3] = mass1[(i + 1) * (ny + 1) + j] - mass1[i * (ny + 1) + j];
				b[1] = 0;
				b[2] = cellsize / (float) ky;
				b[3] = mass1[i * (ny + 1) + j + 1] - mass1[i * (ny + 1) + j];
				c[1] = cellsize / (float) kx;
				c[2] = - cellsize / (float) ky;
				c[3] = mass1[(i + 1) * (ny + 1) + j] - mass1[i * (ny + 1) + j + 1];
				normal[k * 18 + 4 * 3 + 0] = a[2] * b[3] - a[3] * b[2];
				normal[k * 18 + 4 * 3 + 1] = a[3] * b[1] - a[1] * b[3];
				normal[k * 18 + 4 * 3 + 2] = a[1] * b[2] - a[2] * b[1];
				a[0] = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
				b[0] = sqrt(b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
				c[0] = sqrt(c[1] * c[1] + c[2] * c[2] + c[3] * c[3]);
				p = (a[0] + b[0] + c[0]) / 2;
				area[k * 6 + 4] = sqrt(p * (p - a[0]) * (p - b[0]) * (p - c[0]));
				a[1] = cellsize / (float) kx;
				a[2] = 0;
				a[3] = mass1[(i + 1) * (ny + 1) + j + 1] - mass1[i * (ny + 1) + j + 1];
			b[1] = 0;
				b[2] = cellsize / (float) ky;
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
				area[k * 6 + 5] = area[k * 6 + 0];
				area[k * 6 + 2] = cellsize * cellsize / ((float) kx * (float) kz);
				area[k * 6 + 3] = cellsize * cellsize / ((float) kx * (float) kz);
				area[k * 6 + 0] = cellsize * cellsize / ((float) ky * (float) kz);
				area[k * 6 + 1] = cellsize * cellsize / ((float) ky * (float) kz);
				volume[k] = cellsize * cellsize * cellsize / ((float) kx * (float) ky * (float) kz);
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

void DDT_density_velocity(int i, int j, int k)
{
	int p;
	for (p = 0; p < 3; p++) {
		//printf("p = %d\n", p);
		//printf("A_IND(p, i, j, k) = %d\n", A_IND(p, i, j, k));
		//printf("density(i, j, k) = %f\n", density(i, j, k));
		A[A_IND(p, i, j, k)] += density(i, j, k) / dt;
		//printf("A[A_IND(p, i, j, k)] += density(i, j, k) / dt;");
		B[B_IND(p, i, j, k)] += density(i, j, k) * velocity[VEL_IND(p, i, j, k)] / dt;
		//printf("B[B_IND(p, i, j, k)] += density(i, j, k) * velocity[VEL_IND(p, i, j, k)] / dt;");
	}
}

float velocity_on_face(int p, int i, int j, int k, int s)
{
	//printf("velocity_on_face\n");
	if ((ind_cell_multipl[i * ny + j] == -1) ||
		(i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0))
			return 0;
	int a, b;
	//printf("p = %d, i = %d, j = %d, k = %d\n", p, i, j, k);
	switch (s) {
		case 0:
			if ((ind_cell_multipl[(i + 1) * ny + j] == -1) || (i + 1 >= nx))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i + 1, j, k)]) / 2;
		case 1:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
			{
				//printf("p = %d, i = %d, j = %d, k = %d\n", p, i, j, k);
				//printf("velocity[VEL_IND(p, i, j, k)] = %f\nend\n", velocity[VEL_IND(p, i, j, k)]);
				return velocity[VEL_IND(p, i, j, k)];
			}
			else
			{
				//printf("p = %d, i = %d, j = %d, k = %d\n", p, i, j, k);
				//printf("velocity[VEL_IND(p, i, j, k)] = %f\n", velocity[VEL_IND(p, i, j, k)]);
				//printf("VEL_IND(p, i, j, k) = %d\n", VEL_IND(p, i, j, k));
				//printf("VEL_IND(p, i - 1, j, k) = %d\n", VEL_IND(p, i - 1, j, k));
				//printf("PRESS_IND(i, j, k) = %d\n", PRESS_IND(i, j, k));
				//printf("PRESS_IND(i, j + 1, k) = %d\n", PRESS_IND(i, j + 1, k));
				//printf("PRESS_IND(i - 1, j, k) = %d\n", PRESS_IND((i - 1), j, k));
				//printf("PRESS_IND(i - 1, j, k) = %d\n", PRESS_IND(0, 8, 0));
				//printf("(k * n_cells_multipl + ind_cell_multipl[(i - 1) * ny + j]) = %d\n", (k * n_cells_multipl + ind_cell_multipl[(i - 1) * ny + j]));
				//printf("ind_cell_multipl[i * ny + j] = %d\n", ind_cell_multipl[i * ny + j]);
				//printf("ind_cell_multipl[i * ny + j - 1] = %d\n", ind_cell_multipl[i * ny + j - 1]);
				//printf("n_cells_multipl = %d\nind_cell_multipl[(i - 1) * ny + j] = %d\nnx = %d\nny = %d\nnz = %d\n", n_cells_multipl, ind_cell_multipl[(i - 1) * ny + j], nx, ny, nz);
				//for (a = 0; a < nx; a++)
				//	for (b = 0; b < ny; b++)
				//		printf("ind_cell_multipl[%d * ny + %d] = %d\n", a, b, ind_cell_multipl[a * ny + b]);
				//printf("velocity[VEL_IND(p, i - 1, j, k)] = %f\nend\n", velocity[VEL_IND(p, i - 1, j, k)]);
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i - 1, j, k)]) / 2;
			}
		case 2:
			if ((ind_cell_multipl[i * ny + j + 1] == -1) || (j + 1 >= ny))
				return velocity[VEL_IND(p, i, j, k)];
			else
				return (velocity[VEL_IND(p, i, j, k)] + velocity[VEL_IND(p, i, j + 1, k)]) / 2;
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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

float density_on_face(int i, int j, int k, int s)
{
	//printf("density_on_face\n");
	if ((ind_cell_multipl[i * ny + j] == -1) ||
		(i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0))
			return 0;
	switch (s) {
		case 0:
			if ((ind_cell_multipl[(i + 1) * ny + j] == -1) || (i + 1 >= nx))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i + 1, j, k)) / 2;
		case 1:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i - 1, j, k)) / 2;
		case 2:
			if ((ind_cell_multipl[i * ny + j + 1] == -1) || (j + 1 >= ny))
				return density(i, j, k);
			else
				return (density(i, j, k) + density(i, j + 1, k)) / 2;
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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
	int s, p, pr;
	for (pr = 0; pr < 3; pr++) {
		for (s = 0; s < 6; s++) {
			for (p = 0; p < 3; p++) {
				//printf("pr = %d, s = %d, p = %d\n", pr, s, p);
				B[B_IND(p, i, j, k)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
					velocity_on_face(p, i, j, k, s) * normal[NORMAL_IND(p, i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)]);
			}
		}
	}
	return;
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

int A_IND_S_SWITCH(int i, int j, int k, int s)
{
	if ((ind_cell_multipl[i * ny + j] == -1) ||
		(i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0))
			return -1;
	switch (s) {
		case 0:
			if ((ind_cell_multipl[(i + 1) * ny + j] == -1) || (i + 1 >= nx))
				return 0;
			else
				return 1;
		case 1:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
				return 0;
			else
				return 1;
		case 2:
			if ((ind_cell_multipl[i * ny + j + 1] == -1) || (j + 1 >= ny))
				return 0;
			else
				return 1;
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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
	int s, p, pr;
	for (pr = 0; pr < 3; pr++) {
		for (s = 0; s < 6; s++) {
			for (p = 0; p < 3; p++) {
				//printf("1div\n");
				A[A_IND(p, i, j, k)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
					normal[NORMAL_IND(p, i, j, k, s)] / ( 2 * 2 * volume[VOLUME_IND(i, j, k)]);
				//printf("2div\n");
				if (A_IND_S_SWITCH(i, j, k, s) == 1)
				{
					//printf("3div\n");
					//printf("A_IND_S(p, i, j, k, s) = %d\n", A_IND_S(p, i, j, k, s));
					A[A_IND_S(p, i, j, k, s)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
						normal[NORMAL_IND(p, i, j, k, s)] / (2 * 2 * volume[VOLUME_IND(i, j, k)]);
					//printf("3div\n");
				}
				if (A_IND_S_SWITCH(i, j, k, s) == 0)
				{
					//printf("4div\n");
					A[A_IND(p, i, j, k)] -= area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) *
						normal[NORMAL_IND(p, i, j, k, s)] / (2 * 2 * volume[VOLUME_IND(i, j, k)]);
					//printf("4div\n");
				}
				//printf("5div\n");
			}
		}
	}
	return;
}

void DIV_density_velocity_velocity_crank_nikolson(int i, int j, int k)
{
	//printf("DIV(i, j, k, density_velocity_velocity, half_forward_euler);\n");
	DIV(i, j, k, density_velocity_velocity, half_forward_euler);
	//printf("DIV(i, j, k, density_velocity_velocity, half_backward_euler);\n");
	DIV(i, j, k, density_velocity_velocity, half_backward_euler);
	return;
}

void GRAD_pressure_crank_nikolson(int i, int j, int k)
{
	GRAD(i, j, k, pressure, half_forward_euler);
	GRAD(i, j, k, pressure, half_backward_euler);
	return;
}

void GRAD_pressure_half_forward_euler(int i, int j, int k)
{
	if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0)) {
		B[B_IND(0, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure_atmosphere) / (2 * dx);
	} else {
		B[B_IND(0, i, j, k)] -= (pressure[PRESS_IND(i, j, k)] - pressure[PRESS_IND(i - 1, j, k)]) / (2 * dx);
	}
	if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0)) {
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
	if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0)) {
		A[A_IND(0, i, j, k)] -= 1 / (2 * dx);
		B[B_IND(0, i, j, k)] += pressure_atmosphere / (2 * dx);
	} else {
		A[A_IND(0, i, j, k)] -= 1 / (2 * dx);
		A[A_IND(0, i - 1, j, k)] += 1 / (2 * dx);
	}
	if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0)) {
		A[A_IND(1, i, j, k)] -= 1 / (2 * dy);
		B[B_IND(1, i, j, k)] += pressure_atmosphere / (2 * dy);
	} else {
		A[A_IND(1, i, j, k)] -= 1 / (2 * dy);
		A[A_IND(1, i, j - 1, k)] += 1 / (2 * dy);
	}
	if (k - 1 >= 0) {
		A[A_IND(2, i, j, k)] -= 1 / (2 * dz);
		A[A_IND(2, i, j, k - 1)] += 1 / (2 * dz);
	}
	return;
}

void VECT_gravity_force_half_forward_euler(int i, int j, int k)
{
	int p;
	for (p = 0; p < 3; p++) {
		B[B_IND(p, i, j, k)] -= density_snow * g[p] * phase_fraction[PHASE_IND(i, j, k)] / 2 + density_air * g[p] * (1 - phase_fraction[PHASE_IND(i, j, k)]) / 2;
	}
	return;
}

void VECT_gravity_force_half_backward_euler(int i, int j, int k)
{
	int p;
	for (p = 0; p < 3; p++) {
		A[A_IND(3, i, j, k)] += density_snow * g[p] / 2 - density_air * g[p] / 2;
		B[B_IND(p, i, j, k)] -= density_air * g[p] / 2;
	}
	return;
}

void VECT_gravity_force_crank_nikolson(int i, int j, int k)
{
	VECT(i, j, k, gravity_force, half_forward_euler);
	VECT(i, j, k, gravity_force, half_backward_euler);
	return;
}

float strain_rate_on_face(int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
				return 0;
			else
				return (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		case 1:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
				if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
					return 0;
				else
					return 0.5 * (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / dx;
			else
				if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
					return 0.5 * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j - 1, k, s)) / dy;
				else
					return 0.5 * ((velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j - 1, k, s)) / dy + (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / dx);
		case 2:
			if (k - 1 < 0)
				if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
					return 0.5 * velocity_on_face(0, i, j, k, s) / dz;
				else
					return 0.5 * (velocity_on_face(0, i, j, k, s) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / dx);
			else
				if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
					return 0.5 * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j, k - 1, s)) / dz;
				else
					return 0.5 * ((velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i, j, k - 1, s)) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / dx);
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
				return 0;
			else
				return (velocity_on_face(1, i, j, k, s) - velocity_on_face(1, i, j - 1, k, s)) / dx;
		case 5:
			if (k - 1 < 0)
				if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
					return 0.5 * velocity_on_face(1, i, j, k, s) / dz;
				else
					return 0.5 * (velocity_on_face(1, i, j, k, s) / dz + (velocity_on_face(2, i, j, k, s) - velocity_on_face(2, i, j - 1, k, s)) / dy);
			else
				if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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

float shear_rate_on_face(int i, int j, int k, int s)
{
	int m, n;
	float x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate_on_face(i, j, k, s, m, n) * strain_rate_on_face(i, j, k, s, m, n);
		}
	}
	return x;
}

float phase_fraction_on_face(int i, int j, int k, int s)
{
	if ((ind_cell_multipl[i * ny + j] == -1) ||
		(i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0))
			return 1;
	switch (s) {
		case 0:
			if ((ind_cell_multipl[(i + 1) * ny + j] == -1) || (i + 1 >= nx))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i + 1, j, k)]) / 2;
		case 1:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i - 1, j, k)]) / 2;
		case 2:
			if ((ind_cell_multipl[i * ny + j + 1] == -1) || (j + 1 >= ny))
				return phase_fraction[PHASE_IND(i, j, k)];
			else
				return (phase_fraction[PHASE_IND(i, j, k)] + phase_fraction[PHASE_IND(i, j + 1, k)]) / 2;
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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

float effective_viscosity_on_face(int i, int j, int k, int s)
{
	float phase_fraction_on_face_tmp = phase_fraction_on_face(i, j, k, s);
	float shear_rate_on_face_tmp = shear_rate_on_face(i, j, k, s);
	float x = k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= shear_rate_0) {
		x += limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (k_viscosity_snow * pow(shear_rate_on_face_tmp, flow_index - 1) + yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

void DIV_shear_stress_forward_euler(int i, int j, int k)
{
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

float pressure_on_face(int i, int j, int k, int s)
{
	if ((ind_cell_multipl[i * ny + j] == -1) ||
		(i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0))
			return 1;
	switch (s) {
		case 0:
			if ((ind_cell_multipl[(i + 1) * ny + j] == -1) || (i + 1 >= nx))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i + 1, j, k)]) / 2;
		case 1:
			if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i - 1, j, k)]) / 2;
		case 2:
			if ((ind_cell_multipl[i * ny + j + 1] == -1) || (j + 1 >= ny))
				return (pressure[PRESS_IND(i, j, k)] + pressure_atmosphere) / 2;
			else
				return (pressure[PRESS_IND(i, j, k)] + pressure[PRESS_IND(i, j + 1, k)]) / 2;
		case 3:
			if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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
	//printf("DIV(i, j, k, grad_pressure, half_forward_euler);\n");
	DIV(i, j, k, grad_pressure, half_forward_euler);
	//printf("DIV(i, j, k, grad_pressure, half_backward_euler);\n");
	DIV(i, j, k, grad_pressure, half_backward_euler);
	return;
}

void DIV_grad_pressure_half_forward_euler(int i, int j, int k)
{
	int s;
	for (s = 0; s < 6; s++) {
		if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_atmosphere) * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		else
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * (pressure_on_face(i, j, k, s) - pressure_on_face(i - 1, j, k, s)) * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))
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
	int s;
	for (s = 0; s < 6; s++) {
		//printf("x\n");
		if ((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0)) {
			A[A_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			if (A_IND_S_SWITCH(i, j, k, s))
				A[A_IND_S(4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		} else {
			A[A_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			A[A_IND(4, i - 1, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			if (A_IND_S_SWITCH(i - 1, j, k, s) == 1)
				A[A_IND_S(4, i - 1, j, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(0, i, j, k, s)] / dx;
			else
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(0, i, j, k, s)] / dx;
		}
		//printf("y\n");
		if ((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0)) {
			A[A_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			if (A_IND_S_SWITCH(i, j, k, s))
				A[A_IND_S(4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		} else {
			A[A_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			A[A_IND(4, i, j - 1, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			else
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			if (A_IND_S_SWITCH(i, j - 1, k, s) == 1)
				A[A_IND_S(4, i, j - 1, k, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(1, i, j, k, s)] / dy;
			else
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(1, i, j, k, s)] / dy;
		}
		//printf("z\n");
		if (k > 0) {
			A[A_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			//printf("1\n");
			//printf("A_IND(4, i - 1, j, k) = %d\n", A_IND(4, i - 1, j, k));
			A[A_IND(4, i, j, k - 1)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
			//printf("2\n");
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
			{
				A[A_IND_S(4, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				//printf("3\n");
			}
			else
			{
				B[B_IND(4, i, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				//printf("4\n");
			}
			if (A_IND_S_SWITCH(i, j, k - 1, s) == 1)
			{
				A[A_IND_S(4, i, j, k - 1, s)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				//printf("5\n");
			}
			else
			{
				B[B_IND(4, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * pressure_atmosphere * normal[NORMAL_IND(2, i, j, k, s)] / dz;
				//printf("6\n");
			}
		}
	}
	return;
}

void DIV_div_density_velocity_velocity_crank_nikolson(int i, int j, int k)
{
	//printf("DIV(i, j, k, div_density_velocity_velocity, half_forward_euler);\n");
	DIV(i, j, k, div_density_velocity_velocity, half_forward_euler);
	//printf("DIV(i, j, k, div_density_velocity_velocity, half_backward_euler);\n");
	DIV(i, j, k, div_density_velocity_velocity, half_backward_euler);
	return;
}

void DIV_div_density_velocity_velocity_half_forward_euler(int i, int j, int k)
{
	int s;
	for (s = 0; s < 6; s++) {
		if (!((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0)))
			B[B_IND(4, i, j, k)] -= (1 / (2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (velocity_on_face(0, i, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / dx;
		if (!((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0)))
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
	int s;
	for (s = 0; s < 6; s++) {
		if (!((ind_cell_multipl[(i - 1) * ny + j] == -1) || (i - 1 < 0))) {
			A[A_IND(0, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			A[A_IND(0, i - 1, j, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(0, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			else
				A[A_IND(0, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			if (A_IND_S_SWITCH(i - 1, j, k, s) == 1)
				A[A_IND_S(0, i - 1, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
			else
				A[A_IND(0, i - 1, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dx;
		}
		if (!((ind_cell_multipl[i * ny + j - 1] == -1) || (j - 1 < 0))) {
			A[A_IND(1, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			A[A_IND(1, i, j - 1, k)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(1, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			else
				A[A_IND(1, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			if (A_IND_S_SWITCH(i, j - 1, k, s) == 1)
				A[A_IND_S(1, i, j - 1, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
			else
				A[A_IND(1, i, j - 1, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dy;
		}
		if (k == 0) {
			A[A_IND(2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(2, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			else
				A[A_IND(2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
		}
		if (k > 0) {
			A[A_IND(2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			A[A_IND(2, i, j, k - 1)] -= (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			if (A_IND_S_SWITCH(i, j, k, s) == 1)
				A[A_IND_S(2, i, j, k, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			else
				A[A_IND(2, i, j, k)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			if (A_IND_S_SWITCH(i, j, k - 1, s) == 1)
				A[A_IND_S(2, i, j, k - 1, s)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
			else
				A[A_IND(2, i, j, k - 1)] += (1 / (2 * 2 * volume[VOLUME_IND(i, j, k)])) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) * 
					(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
					velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) / dz;
		}
	}
	return;
}

void DDT_density_snow_volume_fraction(int i, int j, int k)
{
	A[A_IND(3, i, j, k)] += density(i, j, k) / dt;
	B[B_IND(3, i, j, k)] += density(i, j, k) * phase_fraction[PHASE_IND(i, j, k)] / dt;
}

void DIV_density_snow_volume_fraction_velocity_crank_nikolson(int i, int j, int k)
{
	//printf("DIV(i, j, k, density_snow_volume_fraction_velocity, half_forward_euler);\n");
	DIV(i, j, k, density_snow_volume_fraction_velocity, half_forward_euler);
	//printf("DIV(i, j, k, density_snow_volume_fraction_velocity, half_backward_euler);\n");
	DIV(i, j, k, density_snow_volume_fraction_velocity, half_backward_euler);
}

void DIV_density_snow_volume_fraction_velocity_half_forward_euler(int i, int j, int k)
{
	int s;
	for (s = 0; s < 6; s++) {
		B[B_IND(3, i, j, k)] -= (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2) * phase_fraction_on_face(i, j, k, s);
	}
	return;
}

void DIV_density_snow_volume_fraction_velocity_half_backward_euler(int i, int j, int k)
{
	int s;
	for (s = 0; s < 6; s++) {
		//printf("DIV_density_snow_volume_fraction_velocity_half_backward_euler 0\n");
		//printf("A[A_IND(3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *\n"
		//	"(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +\n"
		//	 "velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);\n");
		//printf("A[A_IND(3, %d, %d, %d)] += (1 / volume[VOLUME_IND(%d, %d, %d)]) * area[AREA_IND(%d, %d, %d, %d)] * density_on_face(%d, %d, %d, %d) *\n"
		//	"(velocity_on_face(0, %d, %d, %d, %d) * normal[NORMAL_IND(0, %d, %d, %d, %d)] + velocity_on_face(1, %d, %d, %d, %d) * normal[NORMAL_IND(1, %d, %d, %d, %d)] +\n"
		//	 "velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);\n");
		//printf("A[A_IND(3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *\n"
		//	"(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +\n"
		//	 "velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);\n");
		A[A_IND(3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);
		if (A_IND_S_SWITCH(i, j, k, s) == 1) {
			//printf("DIV_density_snow_volume_fraction_velocity_half_backward_euler 1\n");
			//printf("s = %d\n", s);
			//printf("volume[VOLUME_IND(i, j, k)] = %f\n", volume[VOLUME_IND(i, j, k)]);
			//printf("area[AREA_IND(i, j, k, s)] = %f\n", area[AREA_IND(i, j, k, s)]);
			//printf("density_on_face(i, j, k, s) = %f\n", density_on_face(i, j, k, s));
			//printf("velocity_on_face(0, i, j, k, s) = %f\n", velocity_on_face(0, i, j, k, s));
			//printf("normal[NORMAL_IND(0, i, j, k, s)] = %f\n", normal[NORMAL_IND(0, i, j, k, s)]);
			//printf("velocity_on_face(1, i, j, k, s) = %f\n", velocity_on_face(1, i, j, k, s));
			//printf("normal[NORMAL_IND(1, i, j, k, s)] = %f\n", normal[NORMAL_IND(1, i, j, k, s)]);
			//printf("velocity_on_face(2, i, j, k, s) = %f\n", velocity_on_face(2, i, j, k, s));
			//printf("normal[NORMAL_IND(2, i, j, k, s)]) = %f\n", normal[NORMAL_IND(2, i, j, k, s)]);
			//printf("A[A_IND_S(3, i, j, k, s)] = %f\n", A[A_IND_S(3, i, j, k, s)]);
			A[A_IND_S(3, i, j, k, s)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);
			//printf("end\n");
		}
		if (A_IND_S_SWITCH(i, j, k, s) == 0) {
			//printf("DIV_density_snow_volume_fraction_velocity_half_backward_euler 2\n");
			A[A_IND(3, i, j, k)] += (1 / volume[VOLUME_IND(i, j, k)]) * area[AREA_IND(i, j, k, s)] * density_on_face(i, j, k, s) *
				(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] + velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
				 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]) * (1 / 2);
			//printf("end\n");
		}
		//printf("DIV_density_snow_volume_fraction_velocity_half_backward_euler 4\n");
	}
	return;
}

int create_Ab(void)
{
	num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	int system_dimension = n_cells_multipl * nz * num_parameters;
	if ((A = (float *) malloc(system_dimension * system_dimension * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((B = (float *) malloc(system_dimension * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	int i, j, k;
	set_arrays();
	for (k = 0; k < nz; k++) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				if (ind_cell_multipl[i * ny + j] != -1) {
					//printf("create_Ab i = %d,\t\tj = %d,\t\tk = %d\n", i, j, k);
					/*momentum equation*/
					//printf("Adding the momentum equation\n");
					//printf("DDT(i, j, k, density_velocity);\n");
					DDT(i, j, k, density_velocity);
					//printf("DIV(i, j, k, density_velocity_velocity, crank_nikolson);\n");
					DIV(i, j, k, density_velocity_velocity, crank_nikolson);
					//printf("VECT(i, j, k, gravity_force, crank_nikolson);\n");
					VECT(i, j, k, gravity_force, crank_nikolson);
					//printf("GRAD(i, j, k, pressure, crank_nikolson);\n");
					GRAD(i, j, k, pressure, crank_nikolson);
					//printf("DIV(i, j, k, shear_stress, forward_euler);\n");
					DIV(i, j, k, shear_stress, forward_euler);
					/*pressure equation*/
					//printf("Adding the pressure equation\n");
					//printf("DIV(i, j, k, grad_pressure, crank_nikolson);\n");
					DIV(i, j, k, grad_pressure, crank_nikolson);
					//printf("DIV(i, j, k, div_density_velocity_velocity, crank_nikolson);\n");
					DIV(i, j, k, div_density_velocity_velocity, crank_nikolson);
					/*snow volume fraction equation*/
					//printf("Adding the snow volume fraction equation\n");
					//printf("DDT(i, j, k, density_snow_volume_fraction);\n");
					DDT(i, j, k, density_snow_volume_fraction);
					//printf("DIV(i, j, k, density_snow_volume_fraction_velocity, crank_nikolson);\n");
					DIV(i, j, k, density_snow_volume_fraction_velocity, crank_nikolson);
					//printf("End cycle of create_Ab\n");
				}
			}
		}
	}
	//printf("End creating matrix\n");
	//printf("Openning file to write matrix A\n");
	FILE *f = fopen("A.txt","w");
	for (i = 0; i < system_dimension; i++) {
		for (j = 0; j < system_dimension; j++) {
			fprintf(f, "%10lf\t", A[i * system_dimension + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	//printf("Closing file to write matrix A\n");
	return 0;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -a, -s, -p, -v, -k, -i, -l, -S must be declared\n\n");
	printf("-m\tmap\tASCII file of map\n");
	printf("-r\trerion\tASCII file of region\n");
	printf("-H\thight\thight of calculation area (float)\n");
	printf("-D\tdepth\tdepth of snow cover (float)\n");
	printf("-x\tkx\treduction ratio in x direction (int)\n");
	printf("-y\tky\treduction ratio in y direction (int)\n");
	printf("-z\tkz\treduction ratio in z direction (int)\n");
	printf("-a\tsnow density air density\tdensity of snow and air (float)\n");
	printf("-s\tsnow density air density\tdensity of snow and air (float)\n");
	printf("-p\tpressure\tatmosphere pressure\n");
	//printf("-v\tsnow kinematic viscosity air kinematic viscosity\tkinematic viscisity of snow and air (float)\n");
	//printf("-v\tsnow viscosity air viscosity\tviscisity of snow and air (float)\n");
	printf("-v\tair viscosity\tshear viscosity of the air for Newtonian model (float)\n");
	printf("-k\tsnow consistency index\tconsistency index for the constitutive equation of the Herschel-Bulkley model for snow (float)\n");
	printf("-i\tflow index\tflow index for the constitutive equation of the Herschel-Bulkley model for snow (float)\n");
	printf("-l\tyield stress\tyield stress for the constitutive equation of the Herschel-Bulkley model for snow (float)\n");
	printf("-S\tlimiting shear rate\tlimiting shear rate for the effective viscosity to use Herschel-Bulkley model as a generalized Newtonian fluid model (float)\n");
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
		return 1;
	}
	printf("Parameters entered correctly.\n");

	if (read_asc_and_declare_variables() == 1) goto error;
	printf("Files of map and region have been processed\n");
	if (do_interpolation() == 1) goto error;
	printf("The interpolation have been done\n");
	//printf("snow_region[0] = %d\n", snow_region[0]);
	if (print_vtk() == 1) goto error;
	printf("Map was printed to VTK format\n");
	//print_mass();
	dt = 0.01;//we need to set dt!!!
	printf("Creating Ax = b equation\n");
	if (create_Ab() == 1) goto error;
	printf("Ax = b equation well done\n");
	printf("Free massives\n");
	if (free_massives() == 1) goto error;


	return 0;
error:
	printf("Error\n");
	return 1;
}
