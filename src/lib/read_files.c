#include "init_data.h"
#include "read_files.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int read_asc_and_declare_variables(in *I)
{
	printf("Reading maps.\n");
	//printf("!%100s!\n", I->map_name);
	FILE *f = fopen(I->map_name,"r");
	//FILE *f = fopen("maps/map_for_verification_Ab_21.asc","r");
	if (f == NULL) {
		//printf("!%100s!\n", I->map_name);
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

#if AVALANCHE
	f = fopen(I->region_map_name,"r");
	//f = fopen("maps/map_for_verification_Ab_region_2.asc","r");
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
#endif

	f = fopen("map.txt", "r");
#if AVALANCHE
	f1 = fopen("regions_map.txt", "r");
#endif
	int err;
	char str[20];

	if ((err = fscanf(f, "%s %d", str, &I->ncols)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "ncols") != 0) {
		printf("Error file of map: no \"ncols\" tag\n");
		return 1;
	}

	if ((err = fscanf(f, "%s %d", str, &I->nrows)) == EOF) {
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
	
	if ((err = fscanf(f, "%s %lf", str, &I->cellsize)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "cellsize") != 0) {
		printf("Error file of map: no \"cellsize\" tag\n");
		return 1;
	}

	if ((err = fscanf(f, "%s %lf", str, &I->nodata_value)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "NODATA_value") != 0) {
		printf("Error file of map: no \"NODATA_value\" tag\n");
		return 1;
	}

#if AVALANCHE
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

	if (I->cellsize != cellsize1) {
		printf("Cellsize in both maps need to be the same\n");
		return 1;
	}
#endif

	if ((I->mass = (double *) malloc(I->ncols * I->nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->mass, 0.0, I->ncols * I->nrows * sizeof(double));
	if ((I->ind = (int *) malloc(I->ncols * I->nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->ind, 0, I->ncols * I->nrows * sizeof(int));
	if ((I->bl_cond = (int *) malloc((I->ncols - 1) * (I->nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->bl_cond, 0, (I->ncols - 1) * (I->nrows - 1) * sizeof(int));

	for (i = 0; i < I->ncols * I->nrows; i++) {
		if (err = fscanf(f, "%lf", &I->mass[i]) == EOF) {
			printf("Error file of map\n");
			return 1;
		}
		I->ind[i] = -1;
		if (i < (I->ncols - 1) * (I->nrows - 1))
			I->bl_cond[i] = -1;
	}
	fclose(f);
	remove("map.txt");
	
#if AVALANCHE
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

	if ((I->snow_region = (int *) malloc(I->ncols * I->nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->snow_region, 0, I->ncols * I->nrows * sizeof(int));
	double xlucorner = xllcorner - I->cellsize * I->nrows;
	double ylucorner = yllcorner;
	double xlucorner1 = xllcorner1 - cellsize1 * nrows1;
	double ylucorner1 = yllcorner1;
	for (i = 0; i < I->nrows; i++) {
		for (j = 0; j < I->ncols; j++) {
			if ((xlucorner1 - xlucorner >= 0) && (ylucorner -ylucorner1 >= 0) &&
				(j * I->cellsize + ylucorner >= ylucorner1) &&
				(i * I->cellsize + xlucorner >= xlucorner1) &&
				(j * I->cellsize + ylucorner <= ylucorner1 + ncols1 * I->cellsize) &&
				(i * I->cellsize + xlucorner <= xlucorner1 + nrows1 * I->cellsize)) {
					if (mass_tmp[(i - (int) ((xlucorner1 - xlucorner) /\
						I->cellsize)) * ncols1 + j -\
						(int) ((ylucorner1 - ylucorner) /\
						I->cellsize)] == 1)
							I->snow_region[i * I->ncols + j] = 1;
					else if (mass_tmp[(i - (int) ((xlucorner1 - xlucorner) /\
						I->cellsize)) * ncols1 + j -\
						(int) ((ylucorner1 - ylucorner) / I->cellsize)] == 0)
							I->snow_region[i * I->ncols + j] = 0;
					else
						I->snow_region[i * I->ncols + j] = -1;
			} else {
					I->snow_region[i * I->ncols + j] = -1;
			}
		}
	}
	free(mass_tmp);
#endif

	I->n_points = 0;
	I->n_cells = 0;
	for (i = 0; i < I->nrows; i++) {
		for (j = 0; j < I->ncols; j++) {
			if ((i - 1 >= 0) && (j + 1 < I->ncols) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j + 1] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j + 1] != I->nodata_value)) {
					I->ind[i * I->ncols +j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((j + 1 < I->ncols) && (i + 1 < I->nrows) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j + 1] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j + 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((i + 1 < I->nrows) && (j - 1 >= 0) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j - 1] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j - 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j - 1] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j - 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			end_points:
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(I->ind[i * I->ncols + j] != -1) &&
				(I->ind[i * I->ncols + j - 1] != -1) &&
				(I->ind[(i - 1) * I->ncols + j - 1] != -1) &&
				(I->ind[(i - 1) * I->ncols + j] != -1)) {
					I->bl_cond[(i - 1) * (I->ncols - 1) + j - 1] = I->n_cells;
					I->n_cells++;
			}
			if (I->ind[i * I->ncols + j] == -1)
				I->mass[i * I->ncols + j] = I->nodata_value;
		}
	}
	I->nx = (I->nrows - 1) * I->kx;
	I->ny = (I->ncols - 1) * I->ky;
	I->nz = (int) (I->hight / I->cellsize) * I->kz;
	I->dx[0] = I->dx[1] = I->dx[2] = I->cellsize;
	printf("Maps entered correctly.\n");
	return 0;
}
