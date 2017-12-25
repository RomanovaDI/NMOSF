#include "init_data.h"
#include "read_files.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int read_asc(in *I)
{
#if DEBUG
	printf("Reading maps in process %d.\n", I->my_rank);
#endif
	FILE *f = fopen(I->map_name,"r");
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
#if AVALANCHE
	f = fopen(I->region_map_name,"r");
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->mass, 0.0, I->ncols * I->nrows * sizeof(double));
	if ((I->ind = (int *) malloc(I->ncols * I->nrows * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->ind, 0, I->ncols * I->nrows * sizeof(int));
	if ((I->bl_cond = (int *) malloc((I->ncols - 1) * (I->nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
	I->gl_nx = I->nx = (I->nrows - 1) * I->kx;
	I->gl_ny = I->ny = (I->ncols - 1) * I->ky;
	I->gl_nz = I->nz = (int) (I->hight / I->cellsize) * I->kz;
	I->dx[0] = I->dx[1] = I->dx[2] = I->cellsize;
#if DEBUG
	printf("Maps entered correctly in process %d.\n", I->my_rank);
#endif
	return 0;
}

int declare_variables(in *I)
{
#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Setting up global parameters in process %d\n", I->my_rank);
	printf("Process %d, function %s\n", I->my_rank, __func__);
#endif
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&I->gl_nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->gl_ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->gl_nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(I->dx, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->cellsize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->n_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->n_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Process %d PID %d: I->gl_nx = %d, I->gl_ny = %d, I->gl_nz = %d\n", I->my_rank, getpid(), I->gl_nx, I->gl_ny, I->gl_nz);
	printf("Process %d PID %d: I->dx[0] = %lf, I->dx[1] = %lf, I->dx[2] = %lf\n", I->my_rank, getpid(), I->dx[0], I->dx[1], I->dx[2]);
	printf("Process %d PID %d: I->cellsize = %lf, I->ncols = %d, I->nrows = %d\n", I->my_rank, getpid(), I->cellsize, I->ncols, I->nrows);
	printf("Process %d PID %d: I->n_cells = %d, I->n_points = %d\n", I->my_rank, getpid(), I->n_cells, I->n_points);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	return 0;
}

int read_asc_and_declare_variables(in *I)
{
	if (I->my_rank == 0) if (read_asc(I)) return 1;
	if (I->nproc > 1) if (declare_variables(I)) return 1;
	return 0;
}
