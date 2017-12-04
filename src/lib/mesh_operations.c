#include "init_data.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int do_interpolation(in *I) // is running just on 0 process
{
#if DEBUG
	printf("Interpolating mesh.\n");
#endif
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->ind_multipl, 0, ((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(int));

	int count_ind_multipl = 0;
	if ((I->mass1 = (double *) malloc(((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1) * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}

	if ((e = (double *) malloc(I->ncols * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}

	if ((f = (double *) malloc(I->ncols * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
#if DEBUG
		printf("Interpolation along the y axis have been done\n");
#endif
	}
	
	if ((c = (double *) realloc(c, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}

	if ((e = (double *) realloc(e, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}

	if ((f = (double *) realloc(f, I->nrows * sizeof(double))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
#if DEBUG
		printf("Interpolation along the x axes have been done\n");
#endif
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
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
#if DEBUG
	printf("Make boundary function\n");
#endif
	if ((I->ind_boundary_cells = (int *) malloc((I->nx + 2 * I->stencil_size) * (I->ny + 2 * I->stencil_size) * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
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
#if DEBUG
	char file_name[50];
	sprintf(file_name, "tmp/boundary_%d.txt", I->my_rank);
	FILE *f = fopen(file_name, "w");
	for (i = 0; i < I->nx + 2 * I->stencil_size; i++) {
		for (j = 0; j < I->ny + 2 * I->stencil_size; j++) {
			fprintf(f, "%d ", I->ind_boundary_cells[i * (I->ny + 2 * I->stencil_size) + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
#endif
	return 0;
}

int do_decomposition(in *I)
{
#if DEBUG
	printf("Making decomposition in process %d\n", I->my_rank);
#endif
	int i, j;
	if ((I->ind_proc = (int *) malloc(I->gl_nx * I->gl_ny * sizeof(int))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	for (i = 0; i < I->gl_nx * I->gl_ny; i++)
		I->ind_proc[i] = -1;
	if (I->my_rank == 0) {
	//	I->x_regions = (int) round(sqrt((double) I->nproc * ((double) I->nx / (double) I->ny)));
	//	if (I->x_regions < 1) I->x_regions = 1;
	//	I->y_regions = I->nproc / I->x_regions;
		I->nx = I->num_el_in_x_region = (int) round((double) I->gl_nx / (double) I->x_regions);
		I->ny = I->num_el_in_y_region = (int) round((double) I->gl_ny / (double) I->y_regions);
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				I->ind_proc[i * I->gl_ny + j] = (j / I->num_el_in_y_region) * I->y_regions + i / I->num_el_in_x_region;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(I->ind_proc, I->gl_nx * I->gl_ny, MPI_INT, 0, MPI_COMM_WORLD);
	if (I->my_rank == 0) {
		I->gl_ind_cell_multipl = I->ind_cell_multipl;
	} else {
		if ((I->gl_ind_cell_multipl = (int *) malloc(I->gl_nx * I->gl_ny * sizeof(int))) == NULL) {
			printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
			return 1;
		}
	}
	MPI_Bcast(I->gl_ind_cell_multipl, I->gl_nx * I->gl_ny, MPI_INT, 0, MPI_COMM_WORLD);
	if (I->my_rank == 0)
		I->gl_n_cells_multipl = I->n_cells_multipl;
	MPI_Bcast(&I->gl_n_cells_multipl, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->x_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->y_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->num_el_in_x_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->num_el_in_y_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int i_start, j_start;
	if (I->y_regions > 1) {
		if ((I->my_rank + 1) % I->y_regions == 0) {
			I->ny = I->gl_ny - (I->y_regions - 1) * I->num_el_in_y_region + I->stencil_size;
			j_start = (I->y_regions - 1) * I->num_el_in_y_region - I->stencil_size;
		} else if (I->my_rank % I->y_regions == 0) {
			I->ny += I->stencil_size;
			j_start = 0;
		} else {
			I->ny += 2 * I->stencil_size;
			j_start = (I->my_rank % I->y_regions) * I->num_el_in_y_region - I->stencil_size;
		}
	}
	if (I->x_regions > 1) {
		if (I->my_rank >= I->y_regions * (I->x_regions - 1)) {
			I->nx = I->gl_nx - (I->x_regions - 1) * I->num_el_in_x_region + I->stencil_size;
			i_start = (I->x_regions - 1) * I->num_el_in_x_region - I->stencil_size;
		} else if (I->my_rank < I->x_regions) {
			I->nx += I->stencil_size;
			i_start = 0;
		} else {
			I->nx += 2 * I->stencil_size;
			i_start = (I->my_rank / I->y_regions) * I->num_el_in_x_region - I->stencil_size;
		}
	}
	I->ind_start_region_proc[0] = i_start;
	I->ind_start_region_proc[1] = j_start;
	if ((I->coordinates_of_cell = (int *) malloc(I->nx * I->ny * 2)) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	if ((I->ind_cell_multipl = (int *) malloc(I->nx * I->ny * sizeof(int))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	int ind_write;
	I->n_cells_multipl = 0;
	MPI_Status status;
	for (i = i_start; i < i_start + I->nx; i++) {
		for (j = j_start; j < j_start + I->ny; j++) {
			if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
				I->ind_cell_multipl[(i - i_start) * (I->ny) + j - j_start] = I->n_cells_multipl++;
			}
		}
	}
#if DEBUG
	printf("Process %d, func %s\n", I->my_rank, __func__);
#endif
	if ((I->gl_B = (double *) malloc(I->num_parameters * I->gl_nz * I->gl_n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	return 0;
}
