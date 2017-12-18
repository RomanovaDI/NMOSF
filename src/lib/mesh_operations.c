#include "init_data.h"
#include "utils.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

int do_interpolation(in *I) // is running just on 0 process
{
#if DEBUG
	printf("Interpolating mesh in process %d.\n", I->my_rank);
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
		printf("Interpolation along the y axis have been done in process %d\n", I->my_rank);
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
		printf("Interpolation along the x axes have been done in process %d\n", I->my_rank);
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
	printf("Make boundary function in process %d.\n", I->my_rank);
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
//	char file_name[50];
//	sprintf(file_name, "tmp/boundary_%d.txt", I->my_rank);
//	FILE *f = fopen(file_name, "w");
//	for (i = 0; i < I->nx + 2 * I->stencil_size; i++) {
//		for (j = 0; j < I->ny + 2 * I->stencil_size; j++) {
//			fprintf(f, "%d ", I->ind_boundary_cells[i * (I->ny + 2 * I->stencil_size) + j]);
//		}
//		fprintf(f, "\n");
//	}
//	fclose(f);
#endif
	return 0;
}

int do_decomposition(in *I)
{
#if DEBUG
	printf("Making decomposition in process %d\n", I->my_rank);
#endif
	int i, j, k, l;
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
		if (I->gl_nx - (I->x_regions - 1) * I->num_el_in_x_region > I->nx)
			I->max_num_el_in_x_region = I->gl_nx - (I->x_regions - 1) * I->num_el_in_x_region;
		else
			I->max_num_el_in_x_region = I->num_el_in_x_region;
		if (I->gl_ny - (I->y_regions - 1) * I->num_el_in_y_region > I->ny)
			I->max_num_el_in_y_region = I->gl_ny - (I->y_regions - 1) * I->num_el_in_y_region;
		else
			I->max_num_el_in_y_region = I->num_el_in_y_region;
		k = l = 0;
		for (i = 0; i < I->gl_nx; i++) {
			if (i / I->num_el_in_x_region < I->x_regions)
				k = i / I->num_el_in_x_region;
			for (j = 0; j < I->gl_ny; j++) {
				if (j / I->num_el_in_y_region < I->y_regions)
					l = j / I->num_el_in_y_region;
				I->ind_proc[i * I->gl_ny + j] = k * I->y_regions + l;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(I->ind_proc, I->gl_nx * I->gl_ny, MPI_INT, 0, MPI_COMM_WORLD);
	if (I->my_rank == 0) {
		I->gl_ind_cell_multipl = I->ind_cell_multipl;
		I->gl_n_cells_multipl = I->n_cells_multipl;
	} else {
		if ((I->gl_ind_cell_multipl = (int *) malloc(I->gl_nx * I->gl_ny * sizeof(int))) == NULL) {
			printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
			return 1;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(I->gl_ind_cell_multipl, I->gl_nx * I->gl_ny, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->gl_n_cells_multipl, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->x_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->y_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->num_el_in_x_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->num_el_in_y_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->max_num_el_in_x_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->max_num_el_in_y_region, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&I->n_points_multipl, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
	} else {
		j_start = 0;
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
	} else {
		i_start = 0;
	}
	if ((I->ind_start_region_proc = (int *) malloc(I->nproc * 2 * sizeof(int))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	I->ind_start_region_proc[0] = i_start;
	I->ind_start_region_proc[1] = j_start;
	if ((I->ind_cell_multipl = (int *) malloc(I->nx * I->ny * sizeof(int))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	I->n_cells_multipl = 0;
	for (i = i_start; i < i_start + I->nx; i++) {
		for (j = j_start; j < j_start + I->ny; j++) {
			if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
				I->ind_cell_multipl[(i - i_start) * (I->ny) + j - j_start] = I->n_cells_multipl++;
			}
		}
	}
	if ((I->gl_B = (double *) malloc(I->num_parameters * I->gl_nz * I->gl_n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	return 0;
}

int reconstruct_src_1(in *I)
{
#if DEBUG
	printf("Reconstructing domain in process %d.\n", I->my_rank);
#endif
	memset(I->gl_B, 0, I->num_parameters * I->gl_nz * I->gl_n_cells_multipl * sizeof(double));
	MPI_Status status;
	double B_tmp[I->num_parameters];
	int i, j, k, p;
	MPI_Barrier(MPI_COMM_WORLD);
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
					if (I->ind_proc[i * I->gl_ny + j] == I->my_rank) {
						if (I->my_rank == 0)
							for (p = 0; p < I->num_parameters; p++)
								I->gl_B[GL_A_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i - I->ind_start_region_proc[0], j - I->ind_start_region_proc[1], k)];
						else {
							p = 0;
							MPI_Send(&I->B_prev[B_IND(I, p, i - I->ind_start_region_proc[0], j - I->ind_start_region_proc[1], k)], I->num_parameters, MPI_DOUBLE, 0, GL_A_IND(I, p, i, j, k), MPI_COMM_WORLD);
						}
					} else if (I->my_rank == 0) {
						p = 0;
						MPI_Recv(B_tmp, I->num_parameters, MPI_DOUBLE, MPI_ANY_SOURCE, GL_A_IND(I, p, i, j, k), MPI_COMM_WORLD, &status);
						for (p = 0; p < I->num_parameters; p ++)
							I->gl_B[GL_A_IND(I, p, i, j, k)] = B_tmp[p];
					}
				}
			}
		}
	}
#if DEBUG
	printf("Finish reconstructing domain in process %d.\n", I->my_rank);
#endif
	return 0;
}

int reconstruct_src_2(in *I)
{
#if DEBUG
	printf("Reconstructing domain in process %d.\n", I->my_rank);
#endif
	memset(I->gl_B, 0, I->num_parameters * I->gl_nz * I->gl_n_cells_multipl * sizeof(double));
	MPI_Status status;
	int num_el = I->max_num_el_in_x_region * I->max_num_el_in_y_region * I->gl_nz * I->num_parameters;
#if DEBUG
	printf("Process %d: I->nx = %d, I->max_num_el_in_x_region = %d, I->num_el_in_x_region = %d, I->ny = %d, I->max_num_el_in_y_region = %d, I->num_el_in_y_region = %d\n",\
		I->my_rank, I->nx, I->max_num_el_in_x_region, I->num_el_in_x_region, I->ny, I->max_num_el_in_y_region, I->num_el_in_y_region);
	printf("Process %d: num_el = %d, num_el * I->nproc = %d\n", I->my_rank, num_el, num_el * I->nproc);
#endif
	double *B_tmp1;
	if ((B_tmp1 = (double *) malloc(num_el * sizeof(double))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	double *B_tmp2;
	if ((B_tmp2 = (double *) malloc(num_el * I->nproc * sizeof(double))) == NULL) {
		printf("Memory error in func %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	int i, j, k, p, m = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if ((I->ind_cell_multipl[i * I->ny + j] != -1) && (I->ind_proc[(i + I->ind_start_region_proc[0]) * I->gl_ny + j + I->ind_start_region_proc[1]] == I->my_rank)) {
					for (p = 0; p < I->num_parameters; p++) {
						B_tmp1[(k * I->max_num_el_in_x_region * I->max_num_el_in_y_region + \
							(i + I->ind_start_region_proc[0] - (I->my_rank / I->y_regions) * I->max_num_el_in_x_region) * I->max_num_el_in_y_region + \
							j + I->ind_start_region_proc[1] - (I->my_rank % I->y_regions) * I->max_num_el_in_y_region) * I->num_parameters + p] =
								I->B_prev[B_IND(I, p, i, j, k)];
					}
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgather(B_tmp1, num_el, MPI_DOUBLE, B_tmp2, num_el, MPI_DOUBLE, MPI_COMM_WORLD);
	int ind_B_tmp2;
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
					for (p = 0; p < I->num_parameters; p++) {
						ind_B_tmp2 = num_el * I->ind_proc[i * I->gl_ny + j] + \
								(k * I->max_num_el_in_x_region * I->max_num_el_in_y_region +
								(i % I->num_el_in_x_region) * I->max_num_el_in_y_region +
								(j % I->num_el_in_y_region)) * I->num_parameters + p;
					//	printf("Process %d: i = %d, j = %d, k = %d, p = %d, ind_B_tmp2 = %d, I->ind_proc[%d] = %d\n", I->my_rank, i, j, k, p, ind_B_tmp2, i * I->gl_ny + j, I->ind_proc[i * I->gl_ny + j]);
						I->gl_B[GL_A_IND(I, p, i, j, k)] =
							B_tmp2[num_el * I->ind_proc[i * I->gl_ny + j] + \
								(k * I->max_num_el_in_x_region * I->max_num_el_in_y_region + \
								(i % I->num_el_in_x_region) * I->max_num_el_in_y_region + \
								(j % I->num_el_in_y_region)) * I->num_parameters + p];
					}
				}
			}
		}
	}
	free(B_tmp2);
	free(B_tmp1);
#if DEBUG
	printf("Finish reconstructing domain in process %d.\n", I->my_rank);
#endif
	return 0;
}

int reconstruct_src(in *I)
{
//	time_t time1, time2;
	double t;
//	double seconds;
//	time1 = time(NULL);
//	int i;
	t = MPI_Wtime();
//	for (i = 0; i < 100; i++)
//	if (reconstruct_src_1(I)) return 1;
	if (reconstruct_src_2(I)) return 1;
	t = MPI_Wtime() - t;
//	time2 = time(NULL);
//	seconds = difftime(time2, time1);
	printf("Time of comunication between processes in process %d PID %d: %lf\n", I->my_rank, getpid(), t);
	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}
