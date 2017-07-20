#include "init_data.h"
#include "utils.h"
#include "vtk_map_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int print_vtk(in *I, int n)
{
	int nn = n;
	int i = 0, j, k, a;
	while (nn > 0) {
		nn /= 10;
		i++;
	}
	//char file_name[8 + i];
	char file_name[20 + i];
	sprintf(file_name, "result/map%d.vtk", n);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "slope\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

	printf("n_points_multipl = %d\n", I->n_points_multipl);
	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	fprintf(f, "POINTS %d double\n", I->n_points_multipl * ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	for (k = 0; k <= (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (i = 0; i < (I->nrows - 1) * I->kx + 1; i ++) {
			for (j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
				if (I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * I->cellsize / (double) I->kx, j * I->cellsize / (double) I->ky,
						I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] + k * I->cellsize / (double) I->kz);
			}
		}
	}
	
	a = 0;
	fprintf(f, "CELLS %d %d\n", I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize),
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) * 9);
	for (k = 0; k < (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
			for (j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
				if ((I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
					(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
					(I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1) &&
					(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1)) {
						fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8,
							k * I->n_points_multipl + I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j - 1],
							k * I->n_points_multipl + I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j - 1],
							k * I->n_points_multipl + I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j],
							k * I->n_points_multipl + I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j],
							(k + 1) * I->n_points_multipl + I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j - 1],
							(k + 1) * I->n_points_multipl + I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j - 1],
							(k + 1) * I->n_points_multipl + I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j],
							(k + 1) * I->n_points_multipl + I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j]);
				}
			}
		}
	}
	fprintf(f, "CELL_TYPES %d\n", I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize));
	for (i = 0; i < I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize); i++) {
		fprintf(f, "%d\n", 12);
	}

	fprintf(f, "CELL_DATA %d\n", I->n_cells_multipl * I->nz);
	fprintf(f, "VECTORS velocity double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\t%20.10f\t%20.10f\n", I->B_prev[B_IND(I, 0, i, j, k)], I->B_prev[B_IND(I, 1, i, j, k)], I->B_prev[B_IND(I, 2, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\n", I->B_prev[B_IND(I, 4, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS snow_volume_fraction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\n", I->B_prev[B_IND(I, 3, i, j, k)]);
			}
		}
	}
	fclose(f);
	sprintf(file_name, "tmp/B_prev%d.txt", n);
	f = fopen(file_name, "w");
	int p;
	for (p = 0; p < I->num_parameters; p++) {
		fprintf(f, "parameter %d\n", p);
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			fprintf(f, "k = %d\n", k);
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if (I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1)
						fprintf(f, "%20.10lf\t", I->B_prev[B_IND(I, p, i, j, k)]);
				}
				fprintf(f, "\n");
			}
		}
	}
	fclose(f);
	return 0;
}

