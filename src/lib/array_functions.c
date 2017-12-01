#include "init_data.h"
#include "utils.h"
#include "array_functions.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int free_massives(in *I)
{
	free(I->B);
#if AVALANCHE
	free(I->area);
	free(I->normal);
	free(I->volume);
#endif
	free(I->ind_cell_multipl);
	free(I->ind_multipl);
#if AVALANCHE
	free(I->snow_region);
#endif
	free(I->Aelem_csr);
	free(I->Ajptr_csr);
	free(I->Aiptr_csr);
	free(I->B_prev);
	free(I->ind_boundary_cells);
	if (I->my_rank == 0) {
		free_initial_arrays(I);
	}
	if (I->nproc > 1) {
		free(I->gl_ind_cell_multipl);
		free(I->ind_proc);
		free(I->gl_B);
	}
	return 0;
}

int free_initial_arrays(in *I)
{
	free(I->bl_cond);
	free(I->ind);
	free(I->mass);
	free(I->mass1);
	return 0;
}

void print_mass(in *I)
{
	int i;
	printf("mass\n");
	for (i = 0; i < I->ncols * I->nrows; i++)
		printf("%f\n", I->mass[i]);
	printf("mass1\n");
	printf("kx = %d\nky = %d\nkz = %d\n", I->kx, I->ky, I->kz);
	for (i = 0; i < ((I->ncols - 1) * I->ky + 1) * ((I->nrows - 1) * I->kx + 1); i++)
		printf("%f\n", I->mass1[i]);
	return;
}

int set_arrays(in *I)
{
	printf("Set arrays of volume, normales and areas for mesh\n");
	int i, j, k, l, m;
	if ((I->volume = (double *) malloc(I->n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->volume, 0, I->n_cells_multipl * sizeof(double));
	if ((I->normal = (double *) malloc(6 * 3 * I->n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->normal, 0, 6 * 3 * I->n_cells_multipl * sizeof(double));
	if ((I->area = (double *) malloc(6 * I->n_cells_multipl * sizeof(double))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	memset((void *) I->area, 0, 6 * I->n_cells_multipl * sizeof(double));
	double a[4], b[4], c[4], p, d;
	k = 0;
	for (i = 0; i < I->nx; i++) {
		for (j = 0; j < I->ny; j++) {
			if (I->ind_cell_multipl[i * I->ny + j] != -1) {
				a[1] = I->cellsize / (double) I->kx;
				a[2] = 0;
				a[3] = I->mass1[(i + 1) * (I->ny + 1) + j] - I->mass1[i * (I->ny + 1) + j];
				b[1] = 0;
				b[2] = I->cellsize / (double) I->ky;
				b[3] = I->mass1[i * (I->ny + 1) + j + 1] - I->mass1[i * (I->ny + 1) + j];
				c[1] = I->cellsize / (double) I->kx;
				c[2] = - I->cellsize / (double) I->ky;
				c[3] = I->mass1[(i + 1) * (I->ny + 1) + j] - I->mass1[i * (I->ny + 1) + j + 1];
				I->normal[k * 18 + 4 * 3 + 0] = a[2] * b[3] - a[3] * b[2];
				I->normal[k * 18 + 4 * 3 + 1] = a[3] * b[1] - a[1] * b[3];
				I->normal[k * 18 + 4 * 3 + 2] = a[1] * b[2] - a[2] * b[1];
				a[0] = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
				b[0] = sqrt(b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
				c[0] = sqrt(c[1] * c[1] + c[2] * c[2] + c[3] * c[3]);
				p = (a[0] + b[0] + c[0]) / 2;
				I->area[k * 6 + 4] = sqrt(p * (p - a[0]) * (p - b[0]) * (p - c[0]));
				a[1] = I->cellsize / (double) I->kx;
				a[2] = 0;
				a[3] = I->mass1[(i + 1) * (I->ny + 1) + j + 1] - I->mass1[i * (I->ny + 1) + j + 1];
				b[1] = 0;
				b[2] = I->cellsize / (double) I->ky;
				b[3] = I->mass1[(i + 1) * (I->ny + 1) + j + 1] - I->mass1[(i + 1) * (I->ny + 1) + j];
				I->normal[k * 18 + 4 * 3 + 0] += a[2] * b[3] - a[3] * b[2];
				I->normal[k * 18 + 4 * 3 + 1] += a[3] * b[1] - a[1] * b[3];
				I->normal[k * 18 + 4 * 3 + 2] += a[1] * b[2] - a[2] * b[1];
				a[0] = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
				b[0] = sqrt(b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
				p = (a[0] + b[0] + c[0]) / 2;
				I->area[k * 6 + 4] += sqrt(p * (p - a[0]) * (p - b[0]) * (p - c[0]));
				d = sqrt(pow(I->normal[k * 18 + 4 * 3 + 0], 2) + pow(I->normal[k * 18 + 4 * 3 + 1], 2) + pow(I->normal[k * 18 + 4 * 3 + 2], 2));
				I->normal[k * 18 + 4 * 3 + 0] /= d;
				I->normal[k * 18 + 4 * 3 + 1] /= d;
				I->normal[k * 18 + 4 * 3 + 2] /= d;
				I->normal[k * 18 + 5 * 3 + 0] = - I->normal[k * 18 + 4 * 3 + 0];
				I->normal[k * 18 + 5 * 3 + 1] = - I->normal[k * 18 + 4 * 3 + 1];
				I->normal[k * 18 + 5 * 3 + 2] = - I->normal[k * 18 + 4 * 3 + 2];
				I->normal[k * 18 + 2 * 3 + 0] = 0;
				I->normal[k * 18 + 2 * 3 + 1] = 1;
				I->normal[k * 18 + 2 * 3 + 2] = 0;
				I->normal[k * 18 + 3 * 3 + 0] = 0;
				I->normal[k * 18 + 3 * 3 + 1] = - 1;
				I->normal[k * 18 + 3 * 3 + 2] = 0;
				I->normal[k * 18 + 0 * 3 + 0] = 1;
				I->normal[k * 18 + 0 * 3 + 1] = 0;
				I->normal[k * 18 + 0 * 3 + 2] = 0;
				I->normal[k * 18 + 1 * 3 + 0] = - 1;
				I->normal[k * 18 + 1 * 3 + 1] = 0;
				I->normal[k * 18 + 1 * 3 + 2] = 0;
				I->area[k * 6 + 2] = I->cellsize * I->cellsize / ((double) I->kx * (double) I->kz);
				I->area[k * 6 + 3] = I->cellsize * I->cellsize / ((double) I->kx * (double) I->kz);
				I->area[k * 6 + 0] = I->cellsize * I->cellsize / ((double) I->ky * (double) I->kz);
				I->area[k * 6 + 1] = I->cellsize * I->cellsize / ((double) I->ky * (double) I->kz);
				I->area[k * 6 + 5] = I->area[k * 6 + 0];
				I->volume[k] = I->cellsize * I->cellsize * I->cellsize / ((double) I->kx * (double) I->ky * (double) I->kz);
				k++;
			}
		}
	}
	FILE *f = fopen("tmp/volume", "w");
	for (i = 0; i < I->nx; i++) {
		for (j = 0; j < I->ny; j++) {
			if (I->ind_cell_multipl[i * I->ny + j] != -1) {
				fprintf(f, "%20.10lf\t", I->volume[VOLUME_IND(I, i, j, k)]);
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen("tmp/area", "w");
	for (i = 0; i < I->nx; i++) {
		for (j = 0; j < I->ny; j++) {
			if (I->ind_cell_multipl[i * I->ny + j] != -1) {
				fprintf(f, "%20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf\t", I->area[AREA_IND(I, i, j, k, 0)], I->area[AREA_IND(I, i, j, k, 1)],
					I->area[AREA_IND(I, i, j, k, 2)], I->area[AREA_IND(I, i, j, k, 3)], I->area[AREA_IND(I, i, j, k, 4)], I->area[AREA_IND(I, i, j, k, 5)]);
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen("tmp/normal", "w");
	for (i = 0; i < I->n_cells_multipl; i++) {
		for (j = 0; j < 18; j++) {
			fprintf(f, "%20.10lf\t", I->normal[i * 18 + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}

int print_B_prev(in *I)
{
	int i, j, k, p;
	FILE *f;
	char file_name[30];
	for (p = 0; p < I->num_parameters; p++) {
		for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
			sprintf(file_name, "tmp/B_prev_p%d_j%d", p, j);
			f = fopen(file_name, "w");
			for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
				for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
					if (I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) {
						fprintf(f, "%f\t", I->B_prev[B_IND(I, p, i, j, k)]);
					}
					fprintf(f, "\n");
				}
			}
			fclose(f);
		}
	}
}

