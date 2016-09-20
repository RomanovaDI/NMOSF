#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int read_asc_and_declare_variables(void);
int do_interpolation(void);
int print_vtk(void);
int free_massives(void);
void display_usage(void);
int create_A(void);
int set_arrays(void);
int annihilate_array(void *a, int size_bites);

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
int n_points_multipl;
int *ind_multipl;
float interpolation, interpolation_poli;
float *mass1;
float nodata_value;
int *snow_region;
float depth;
float density_snow, density_air;
float *phase_fraction;
float *density;
float *velocity;
float *pressure;
float pressure_atmosphere;
float g;

int read_asc_and_declare_variables(void)
{
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

	int *snow_region = (int *) malloc(ncols * nrows * sizeof(int));
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
	return 0;
}

int annihilate_array(void *a, int size_bites)
{
	char *b;
	b = (char *) a;
	int i;
	for (i = 0; i < size_bites; i++)
		b[i] = '0';
}

int do_interpolation(void)
{
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

	if ((ind_multipl = (int *) malloc(ncols * kx * nrows * ky * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	int count_ind_multipl = 0;
	if ((mass1 = (float *) malloc(ncols * kx * nrows * ky * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	for (i = 0; i < ncols * nrows * kx * ky; i++)
		mass1[i] = nodata_value;
	for (i = 0; i < nrows; i++)
		for (j = 0; j < ncols; j++)
			mass1[i * ncols * kx * ky + j * ky] = mass[i * ncols + j];

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
				ind_finish = j + 1;
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
					mass1[i * kx * ncols * ky + (j - 1) * ky + k] = a + b * ((float) k * cellsize / (float) ky - cellsize) +
						c[j] * pow((float) k * cellsize / (float) ky - cellsize, 2) / 2 +
						d * pow((float) k * cellsize / (float) ky - cellsize, 3) / 6;
				}
			}
		}
		annihilate_array((void *) c, ncols * sizeof(float));
		annihilate_array((void *) e, ncols * sizeof(float));
		annihilate_array((void *) f, ncols * sizeof(float));
	}
	
	printf("After interpolation along y axis\n");
	for (i = 0; i < ncols * nrows * kx * ky; i++)
		printf("%f\n", mass1[i]);
	printf("After interpolation along both of axis\n");

	free(c);
	free(e);
	free(f);

	printf("Interpolation along the y axis have been done\n");
	
	if ((c = (float *) malloc(nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (float *) malloc(nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (float *) malloc(nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	for (j = 0; j < ncols; j++) {
		for (l = 0; l < ky; l++) {
			flag = 0;
			ind_start = 0;
			ind_finish = nrows;
			for (i = 0; i < nrows; i++) {
				if ((mass[i * ncols + j] != nodata_value) && (flag == 0)) {
					ind_start = i;
					flag = 1;
				}
				if ((mass[i * ncols + j] == nodata_value) && (flag == 1)) {
					ind_finish = i + 1;
					flag = 2;
				}
			}
			e[ind_start + 2] = f[ind_start + 2] = 0;
			for (i = ind_start + 3; i < ind_finish - 1; i++) {
				e[i] = - 1. / (4 * e[i - 1] + 1);
				f[i] = (6 * (mass1[(i + 1) * ncols * ky * kx + j * ky + l] - 2 * mass1[i * ncols * ky * kx + j * ky + l] + mass1[(i - 1) * ncols * ky * kx + j * ky + l]) /
					(cellsize * cellsize) - 4 * f[i - 1]) / (4 * e[i - 1] + 1);
			}
			i = ind_finish - 2;
			c[i] = (6 * (mass1[(i + 1) * ncols * ky * kx + j * ky + l] - mass1[i * ncols * ky * kx + j * ky + l] + mass1[(i - 1) * ncols * kx * ky + j * ky + l]) /
				(cellsize * cellsize) - 4 * f[i]) / (1 + 4 * e[i]);
			for (i = ind_finish - 3; i > ind_start + 1; i--) {
				c[i] = e[i + 1] * c[i + 1] + f[i + 1];
			}
			c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
			for (i = ind_start + 1; i < ind_finish; i++) {
				a = mass1[i * ncols * ky * kx + j * ky + l];
				d = (c[i] - c[i - 1]) / cellsize;
				b = (mass1[i * ncols * kx * ky + j * ky + l] - mass1[(i - 1) * ncols * kx * ky + j * ky + l]) /
					cellsize + cellsize * (2 * c[i] + c[i - 1]) / 6;
				if ((mass1[(i - 1) * ncols * kx * ky + j * ky + l] != nodata_value) && (mass1[i * ncols * kx * ky + j * ky + l] != nodata_value)) {
						for (k = 0; k < kx; k++) {
							mass1[(i - 1) * kx * ncols * ky + k * ncols * ky + j * ky + l] = a + b * ((float) k * cellsize / (float) kx - cellsize) +
								c[i] * pow((float) k * cellsize / (float) kx - cellsize, 2) / 2 +
								d * pow((float) k * cellsize / (float) kx - cellsize, 3) / 6;
						}
				}
			}
			annihilate_array((void *) c, nrows * sizeof(float));
			annihilate_array((void *) e, nrows * sizeof(float));
			annihilate_array((void *) f, nrows * sizeof(float));
		}
	}

	free(c);
	free(e);
	free(f);

	printf("Interpolation along the x axes have been done\n");

	count_ind_multipl = 0;
	for (i = 0; i < nrows * kx; i++) {
		for (j = 0; j < ncols * ky; j++) {
			if (mass1[i * ncols * ky + j] != nodata_value) {
				ind_multipl[i * ncols * ky + j] = count_ind_multipl;
				count_ind_multipl++;
			} else {
				ind_multipl[i * ncols * ky + j] = -1;
			}
		}
	}

	printf("%d - n_points_multipl\n%d - count_ind_multipl\n", n_points_multipl, count_ind_multipl);

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
		for (i = 0; i < nrows * kx; i ++) {
			for (j = 0; j < ncols * ky; j++) {
				if (ind_multipl[i * ncols * ky + j] != -1)
					fprintf(f, "%f %f %f\n", j * cellsize / (float) kx, i * cellsize / (float) ky,
						mass1[i * ncols * kx + j] + k * cellsize / (float) kz);
			}
		}
	}
	
	a = 0;
	fprintf(f, "CELLS %d %d\n", n_cells * kx * ky * kz * (int) (hight / cellsize) * 5,
		n_cells * kx * ky * kz * (int) (hight / cellsize) * 5 * 5);
	for (k = 0; k < (int) (hight / cellsize) * kz; k++) {
		for (i = 1; i < nrows * ky; i++) {
			for (j = 1; j < ncols * kx; j++) {
				if ((ind_multipl[i * ncols * kx + j] != -1) &&
					(ind_multipl[(i - 1) * ncols * kx + j] != -1) &&
					(ind_multipl[i * ncols * kx + j - 1] != -1) &&
					(ind_multipl[(i - 1) * ncols * kx + j - 1] != -1)) {
//						fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8,
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1],
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1]);
//						fprintf(f, "%d %d %d %d %d\n", 4,
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1],
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1]);
//						fprintf(f, "%d %d %d %d %d\n", 4,
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j],
//							(k + 1) * n_points_multipl + ind_multipl[i * ncols * kx + j]);
//						fprintf(f, "%d %d %d %d %d\n", 4,
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[i * ncols * kx + j]);
//						fprintf(f, "%d %d %d %d %d\n", 4,
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							(k + 1) * n_points_multipl + ind_multipl[i * ncols * kx + j]);
//						fprintf(f, "%d %d %d %d %d\n", 4,
//							k * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j],
//							k * n_points_multipl + ind_multipl[i * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * kx + j - 1],
//							(k + 1) * n_points_multipl + ind_multipl[i * ncols * kx + j]);
						a++;
				}
			}
		}
	}
	fprintf(f, "CELL_TYPES %d\n", n_cells * kx * ky * kz * (int) (hight / cellsize) * 5);
	for (i = 0; i < n_cells * kx * ky * kz * (int) (hight / cellsize) * 5; i++) {
		fprintf(f, "%d\n", 10);
//		fprintf(f, "%d\n", 12);
	}
	return 0;
}

int free_massives(void)
{
	free(mass);
	free(bl_cond);
	free(ind);
	free(ind_multipl);
	free(mass1);
	return 0;
}

int create_A(void)
{}

void print_mass(void)
{
	int i;
	printf("mass\n");
	for (i = 0; i < ncols * nrows; i++)
		printf("%f\n", mass[i]);
	printf("mass1\n");
	printf("kx = %d\nky = %d\nkz = %d\n", kx, ky, kz);
	for (i = 0; i < ncols * ky * nrows * kx; i++)
		printf("%f\n", mass1[i]);
	return;
}

int set_arrays(void)
{
	int i, j, k, l, m, n;
	if ((phase_fraction = (float *) malloc(n_cells * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((density = (float *) malloc(n_cells * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((velocity = (float *) malloc(n_cells * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((pressure = (float *) malloc(n_cells * nz * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	n = 0;
	for (k = 0; k < nz; n++) {
		for (i = 0; i < nrows - 1; i++) {
			for (l = 0; l < kx; l++) {
				for (j = 0; j < ncols - 1; j++) {
					for (m = 0; m < ky; m++) {
						if ((bl_cond[i * (ncols - 1) + j] != -1) &&
							(snow_region[i * ncols + j] == 1)) {
								if (((float) (k + 1) * cellsize > depth) && (depth > (float) k * cellsize)) {
									phase_fraction[n] = (depth - (float) k * cellsize) / (float) cellsize;
									density[n] = phase_fraction[n] * density_snow + (1 - phase_fraction[n]) * density_air;
								}
								else if (depth > (float) (k + 1) * cellsize) {
									phase_fraction[n] = 1;
									density[n] = density_snow;
								}
								else {
									phase_fraction[n] = 0;
									density[n] = density_air;
								}
						}
						if (bl_cond[i * (ncols - 1) + j] != -1) {
							velocity[n] = 0;
							n++;
						}
					}
				}
			}
		}
	}
	return 0;
}

void display_usage(void)
{
	printf("Options -m, -r, -H, -D, -x, -y, -z, -s, -a, -p must be declared\n\n");
	printf("-m\tmap\tASCII file of map\n");
	printf("-r\trerion\tASCII file of region\n");
	printf("-H\thight\thight of calculation area (float)\n");
	printf("-D\tdepth\tdepth of snow cover (float)\n");
	printf("-x\tkx\treduction ratio in x direction (int)\n");
	printf("-y\tky\treduction ratio in y direction (int)\n");
	printf("-z\tkz\treduction ratio in z direction (int)\n");
	printf("-s\tsnow density\tdensity of snow (float)\n");
	printf("-a\tair density\tdensity of air (float)\n");
	printf("-p\tpressure\tatmosphere pressure\n");
	printf("-h\tdisplay usage\n");
}

int main(int argc, char **argv)
{
	int opt = 0;
	g = 9,81;
	static const char *optString = "m:r:H:D:x:y:z:s:a:p:h?";
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
			case 's':
				density_snow = atof(optarg);
				break;
			case 'a':
				density_air = atof(optarg);
				break;
			case 'p':
				pressure_atmosphere = atof(optarg);
				break;
			case 'h':
			case '?':
				display_usage();
				return 1;
		}
	}
	if (argc != 21) {
		printf("Not enouth arguments\n");
		return 1;
	}

	if (read_asc_and_declare_variables() == 1) goto error;
	printf("Files of map and region have been processed\n");
	if (do_interpolation() == 1) goto error;
	printf("The interpolation have been done\n");
	if (print_vtk() == 1) goto error;
	printf("Map was printed to VTK format\n");
	print_mass();
	if (free_massives() == 1) goto error;


	return 0;
error:
	printf("Error\n");
	return 1;
}
