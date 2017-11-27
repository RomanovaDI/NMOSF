#include "init_data.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int do_interpolation(in *I)
{
	//printf("Interpolating mesh.\n");
	int i, j, k, l = 0;
	double a, b, d;
	double *c;
	double *e;
	double *f;
	
	I->n_points_multipl = 0;
	for (i = 0; i < I->nrows - 1; i++) {
		for (j = 0; j < I->ncols - 1; j++) {
			if (I->bl_cond[i * (I->ncols - 1) + j] != -1) {
				if ((j - 1 < 0) || (I->bl_cond[i * (I->ncols - 1) + j - 1] == -1)) {
					I->n_points_multipl += I->kx + 1;
					if (((i - 1 >= 0) && (I->bl_cond[(i - 1) * (I->ncols - 1) + j] != -1)) ||
						((i - 1 >= 0) && (j - 1 >= 0) && (I->bl_cond[(i - 1) * (I->ncols - 1) + j - 1] != -1))) {
							I->n_points_multipl -= 1;
					}
				}
				if ((i - 1 < 0) || (I->bl_cond[(i - 1) * (I->ncols - 1) + j] == -1)) {
					I->n_points_multipl += I->ky;
					if ((i - 1 >= 0) && (j + 1 < I->ncols - 1) && (I->bl_cond[(i - 1) * (I->ncols - 1) + j + 1] != -1)) {
						I->n_points_multipl -= 1;
					}
				}
				I->n_points_multipl += I->kx * I->ky;
			}
		}
	}

	if ((I->ind_multipl = (int *) malloc(((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->ind_multipl, 0, ((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(int));

	int count_ind_multipl = 0;
	if ((I->mass1 = (double *) malloc(((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->mass1, 0, ((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(double));
	for (i = 0; i < ((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1); i++)
		I->mass1[i] = I->nodata_value;
	for (i = 0; i < I->nrows; i++)
		for (j = 0; j < I->ncols; j++)
			I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky] = I->mass[i * I->ncols + j];

	//printf("mass1 after assignment mass\n");
	//for (i = 0; i < ((ncols - 1) * ky + 1) * ((nrows - 1) * kx + 1); i++)
	//	printf("%f\n", mass1[i]);

	if ((c = (double *) malloc(I->ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (double *) malloc(I->ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (double *) malloc(I->ncols * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	memset((void *) c, 0, I->ncols * sizeof(double));
	memset((void *) e, 0, I->ncols * sizeof(double));
	memset((void *) f, 0, I->ncols * sizeof(double));

	int ind_start, ind_finish, flag;

	if ((I->ncols > 2) && (I->ky > 1)) {
		for (i = 0; i < I->nrows; i++) {
			flag = 0;
			ind_start = 0;
			ind_finish = I->ncols;
			for (j = 0; j < I->ncols; j++) {
				if ((I->ind[i * I->ncols + j] != -1) && (flag == 0)) {
					ind_start = j;
					flag = 1;
				}
				if ((I->ind[i * I->ncols + j] == -1) && (flag == 1)) {
					ind_finish = j;
					flag = 2;
				}
			}
			e[ind_start + 2] = f[ind_start + 2] = 0;
			for (j = ind_start + 3; j < ind_finish - 1; j++) {
				e[j] = - 1. / (4 * e[j - 1] + 1);
				f[j] = (6 * (I->mass[i * I->ncols + j + 1] - 2 * I->mass[i * I->ncols + j] + I->mass[i * I->ncols + j - 1]) / (I->cellsize * I->cellsize) - 4 * f[j - 1]) / (4 * e[j - 1] + 1);
			}
			j = ind_finish - 2;
			c[j] = (6 * (I->mass[i * I->ncols + j + 1] - 2 * I->mass[i * I->ncols + j] + I->mass[i * I->ncols + j - 1]) / (I->cellsize * I->cellsize) - 4 * f[j]) / (1 + 4 * e[j]);
			for (j = ind_finish - 3; j > ind_start + 1; j--) {
				c[j] = e[j + 1] * c[j + 1] + f[j + 1];
			}
			c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
			for (j = ind_start + 1; j < ind_finish; j++) {
				a = I->mass[i * I->ncols + j];
				d = (c[j] - c[j - 1]) / I->cellsize;
				b = (I->mass[i * I->ncols + j] - I->mass[i * I->ncols + j - 1]) / I->cellsize + I->cellsize * (2 * c[j] + c[j - 1]) / 6;
				if ((I->ind[i * I->ncols + j - 1] != -1) && (I->ind[i * I->ncols + j] != -1)) {
					for (k = 0; k < I->ky; k++) {
						I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + (j - 1) * I->ky + k] = a + b * (k * I->cellsize / (double) I->ky - I->cellsize) +
							c[j] * pow(k * I->cellsize / (double) I->ky - I->cellsize, 2) / 2 +
							d * pow(k * I->cellsize / (double) I->ky - I->cellsize, 3) / 6;
					}
				}
			}
			memset((void *) c, 0, I->ncols * sizeof(double));
			memset((void *) e, 0, I->ncols * sizeof(double));
			memset((void *) f, 0, I->ncols * sizeof(double));
		}

		printf("Interpolation along the y axis have been done\n");
	}
	
	if ((c = (double *) realloc(c, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((e = (double *) realloc(e, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	if ((f = (double *) realloc(f, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	memset((void *) c, 0, I->nrows * sizeof(double));
	memset((void *) e, 0, I->nrows * sizeof(double));
	memset((void *) f, 0, I->nrows * sizeof(double));

	if ((I->nrows > 2) && (I->kx > 1)) {
		for (j = 0; j < I->ncols; j++) {
			for (l = 0; l < I->ky; l++) {
				flag = 0;
				ind_start = 0;
				ind_finish = I->nrows;
				for (i = 0; i < I->nrows; i++) {
					if (((I->ind[i * I->ncols + j] != -1) && (flag == 0) && (l == 0)) ||
						(((I->ind[i * I->ncols + j] != -1) && ((j + 1 == I->ncols) || (I->ind[i * I->ncols + j + 1] != -1))) &&
						 	(flag == 0) && (l != 0))) {
								ind_start = i;
								flag = 1;
					}
					if (((I->ind[i * I->ncols + j] == -1) && (flag == 1) && (l == 0)) ||
						(((I->ind[i * I->ncols + j] == -1) || ((j + 1 == I->ncols) || (I->ind[i * I->ncols + j + 1] == -1))) &&
						 	(flag == 1) && (l != 0))) {
								ind_finish = i;
								flag = 2;
					}
				}
				e[ind_start + 2] = f[ind_start + 2] = 0;
				for (i = ind_start + 3; i < ind_finish - 1; i++) {
					e[i] = - 1. / (4 * e[i - 1] + 1);
					f[i] = (6 * (I->mass1[(i + 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
						2 * I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] +
						I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
						(I->cellsize * I->cellsize) - 4 * f[i - 1]) / (4 * e[i - 1] + 1);
				}
				i = ind_finish - 2;
				c[i] = (6 * (I->mass1[(i + 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
					2 * I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] +
					I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
					(I->cellsize * I->cellsize) - 4 * f[i]) / (1 + 4 * e[i]);
				for (i = ind_finish - 3; i > ind_start + 1; i--) {
					c[i] = e[i + 1] * c[i + 1] + f[i + 1];
				}
				c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
				for (i = ind_start + 1; i < ind_finish; i++) {
					a = I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l];
					d = (c[i] - c[i - 1]) / I->cellsize;
					b = (I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
						I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
						I->cellsize + I->cellsize * (2 * c[i] + c[i - 1]) / 6;
					//printf("j = %d, l = %d, i = %d, a = %f, b = %f, c = %f, d = %f\n", j, l, i, a, b, c[i], d);
					if ((I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] != I->nodata_value) &&
						(I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] != I->nodata_value)) {
							for (k = 0; k < I->kx; k++) {
								I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + k * ((I->ncols - 1) * I->ky + 1) + j * I->ky + l] =
									a + b * (k * I->cellsize / (double) I->kx - I->cellsize) +
									c[i] * pow(k * I->cellsize / (double) I->kx - I->cellsize, 2) / 2 +
									d * pow(k * I->cellsize / (double) I->kx - I->cellsize, 3) / 6;
							}
					}
				}
				memset((void *) c, 0, I->nrows * sizeof(double));
				memset((void *) e, 0, I->nrows * sizeof(double));
				memset((void *) f, 0, I->nrows * sizeof(double));
				if (j == I->ncols - 1) break;
			}
		}

		printf("Interpolation along the x axes have been done\n");
	}

	free(c);
	free(e);
	free(f);

	count_ind_multipl = 0;
	for (i = 0; i < (I->nrows - 1) * I->kx + 1; i++) {
		for (j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
			if (I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] != I->nodata_value) {
				I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] = count_ind_multipl;
				count_ind_multipl++;
			} else {
				I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] = -1;
			}
		}
	}

	if ((I->ind_cell_multipl = (int *) malloc((I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
//	memset((void *) I->ind_cell_multipl, 0, (I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky * sizeof(int));

	for (i = 0; i < (I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky; i++) {
		I->ind_cell_multipl[i] = -1;
	}

	I->n_cells_multipl = 0;
	for (i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
		for (j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
			if ((I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
				(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
				(I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1) &&
				(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1)) {
					I->ind_cell_multipl[(i - 1) * (I->ncols - 1) * I->ky + j - 1] = I->n_cells_multipl;
					(I->n_cells_multipl)++;
			}
		}
	}

	I->dx[0] /= (double) I->kx;
	I->dx[1] /= (double) I->ky;
	I->dx[2] /= (double) I->kz;
	return 0;
}

int make_boundary(in *I)
{
	printf("Make boundary function\n");
	if ((I->ind_boundary_cells = (int *) malloc((I->nx + 2 * I->stencil_size) * (I->ny + 2 * I->stencil_size) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->ind_boundary_cells, -1, (I->nx + 2 * I->stencil_size) * (I->ny + 2 * I->stencil_size) * sizeof(int));
	int i, j, k, l;
	for (i = I->stencil_size; i < I->stencil_size + I->nx; i++) {
		for (j = I->stencil_size; j < I->stencil_size + I->ny; j++) {
			if (I->ind_cell_multipl[(i - I->stencil_size) * I->ny + j - I->stencil_size] != -1) {
				for (k = -I->stencil_size; k <= I->stencil_size; k++) {
					for (l = -I->stencil_size; l<= I->stencil_size; l++) {
						I->ind_boundary_cells[(i + k) * (I->ny + 2 * I->stencil_size) + j + l] = 1;
					}
				}
			}
		}
	}
	k = 0;
	for (i = 0; i < I->nx + 2 * I->stencil_size; i++) {
		for (j = 0; j < I->ny + 2 * I->stencil_size; j++) {
			if (I->ind_boundary_cells[i * (I->ny + 2 * I->stencil_size) + j] != -1) {
				I->ind_boundary_cells[i * (I->ny + 2 * I->stencil_size) + j] = k;
				k++;
			}
		}
	}
	I->n_boundary_cells = k;
	FILE *f = fopen("tmp/boundary.txt", "w");
	for (i = 0; i < I->nx + 2 * I->stencil_size; i++) {
		for (j = 0; j < I->ny + 2 * I->stencil_size; j++) {
			fprintf(f, "%d ", I->ind_boundary_cells[i * (I->ny + 2 * I->stencil_size) + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}

int do_decomposition(in *I)
{
	int *boundary_marker;
	if (I->my_rank == 0) {
		printf("Making decomposition\n");
		int i, j, num_el_in_x_region, num_el_in_y_region;
		I->x_regions = (int) round(sqrt((double) I->nproc * ((double) I->nx / (double) I->ny)));
		if (I->x_regions < 1) I->x_regions = 1;
		I->y_regions = I->nproc / I->x_regions;
		I->nx = num_el_in_x_region = (int) round((double) I->gl_nx / (double) I->x_regions);
		I->ny = num_el_in_y_region = (int) round((double) I->gl_ny / (double) I->y_regions);
		if ((I->ind_proc = (int *) malloc(2 * I->gl_nx * I->gl_ny * sizeof(int))) == NULL) {
			printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
			return 1;
		}
		if ((boundary_marker = (int *) malloc(I->nproc * 4 * sizeof(int))) == NULL) {
			printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
			return 1;
		}
		memset(boundary_marker, 0, I->nproc * 4 * sizeof(int));
		int ind_cell_proc[I->nproc];
		memset(ind_cell_proc, 0, I->nproc * sizeof(int));
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				I->ind_proc[i * I->gl_ny + j] = (j / num_el_in_y_region) * I->y_regions + i / num_el_in_x_region;
				if (I->ind_cell_multipl[i * I->gl_ny + j] != -1) {
					I->ind_proc[2 * I->gl_nx * I->gl_ny + i * I->gl_ny + j] = ind_cell_proc[(j / num_el_in_y_region) * I->y_regions + i / num_el_in_x_region]++;
				} else {
					I->ind_proc[2 * I->gl_nx * I->gl_ny + i * I->gl_ny + j] = -1;
				}
				if (i == I->gl_nx - 1) boundary_marker[I->ind_proc[i * I->gl_ny + j]] += 1;
				if (i == 0) boundary_marker[I->ind_proc[i * I->gl_ny + j] + 1] += 1;
				if (j == I->gl_ny - 1) boundary_marker[I->ind_proc[i * I->gl_ny + j] + 2] += 1;
				if (j == 0) boundary_marker[I->ind_proc[i * I->gl_ny + j] + 3] += 1;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(I->ind_proc, 2 * I->gl_nx * I->gl_ny, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(boundary_marker, I->nproc * 4, MPI_INT, 0, MPI_COMM_WORLD);
	for (i = 0; i < 4; i++)
		I->boundary_marker[i] = boundary_marker[I->my_rank + i];
	if (I->my_rank == 0) {
		I->gl_ind_cell_multipl = I->ind_cell_multipl;
	}
	MPI_Bcast(I->x_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(I->y_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(I->nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(I->ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(I->nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (I->my_rank + 1 % I->y_regions == 0) {
		I->ny = I->gl_ny - (I->y_regions - 1) * I->ny;
	}
	if (I->my_rank >= I->y_regions * (I->x_regions - 1)) {
		I->nx = I->gl_nx - (I->x_regions - 1) * I->nx;
	}
	return 0;
}
