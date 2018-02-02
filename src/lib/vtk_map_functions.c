#include "init_data.h"
#include "utils.h"
#include "vtk_map_functions.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int print_vtk_avalanche(in *I, int n)
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
//	printf("n_points_multipl = %d\n", I->n_points_multipl);
//	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
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

int print_vtk_termogas(in *I, int n)
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
//	printf("n_points_multipl = %d\n", I->n_points_multipl);
//	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
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
//	fprintf(f, "COLOR_SCALARS concentration 4\n");
//	for (k = 0; k < I->nz; k++) {
//		for (i = 0; i < I->nx; i++) {
//			for (j = 0; j < I->ny; j++) {
//				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
//					for (p = 0; p < 4; p++)
//						fprintf(f, "%20.10f\t", I->B_prev[B_IND(I, p, i, j, k)]);
//					fprintf(f, "\n");
//				}
//			}
//		}
//	}
	int p;
	for (p = 0; p < 4; p++) {
		fprintf(f, "SCALARS concentration_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (k = 0; k < I->nz; k++) {
			for (i = 0; i < I->nx; i++) {
				for (j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						fprintf(f, "%20.20f\n", I->B_prev[B_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", I->B_prev[B_IND(I, 4, i, j, k)]);
			}
		}
	}
//	fprintf(f, "COLOR_SCALARS saturation 3\n");
//	for (k = 0; k < I->nz; k++) {
//		for (i = 0; i < I->nx; i++) {
//			for (j = 0; j < I->ny; j++) {
//				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
//					for (p = 5; p < 8; p++)
//						fprintf(f, "%20.10f\t", I->B_prev[B_IND(I, p, i, j, k)]);
//					fprintf(f, "\n");
//				}
//			}
//		}
//	}
	for (p = 5; p < 8; p++) {
		fprintf(f, "SCALARS saturation_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (k = 0; k < I->nz; k++) {
			for (i = 0; i < I->nx; i++) {
				for (j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						fprintf(f, "%20.20f\n", I->B_prev[B_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS temperature_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", I->B_prev[B_IND(I, 8, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS temperature_environment double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
//	fprintf(f, "SCALARS temperature_flow double 1\n");
//	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", I->B_prev[B_IND(I, 9, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS viscosity_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS density_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS density_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS density_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS rate_of_reaction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", rate_of_reaction(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	double tmp;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					tmp = (concentration(I, 1, i, j, k) * saturation(I, 2, i, j, k) *
						(pressure(I, i, j, k) * I->molar_weight[1] / (I->R * temperature_flow(I, i, j, k))) / I->molar_weight[1]);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_oxigen double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					tmp = (saturation(I, 1, i, j, k) * density_t(I, 1, i, j, k) / I->molar_weight[4]);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS chemical_reaction_heat_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", chemical_reaction_heat_flow(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_water double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_oil double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 1, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_gas double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fclose(f);
	return 0;
}

int print_vtk_termogas_parallel(in *I, int n)
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
#if DEBUG
	printf("n_points_multipl = %d\n", I->n_points_multipl);
	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
#endif
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
	fprintf(f, "CELL_DATA %d\n", I->gl_n_cells_multipl * I->gl_nz);
//	fprintf(f, "COLOR_SCALARS concentration 4\n");
//	for (k = 0; k < I->gl_nz; k++) {
//		for (i = 0; i < I->gl_nx; i++) {
//			for (j = 0; j < I->gl_ny; j++) {
//				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
//					for (p = 0; p < 4; p++)
//						fprintf(f, "%20.10f\t", I->gl_B[GL_A_IND(I, p, i, j, k)]);
//					fprintf(f, "\n");
//				}
//			}
//		}
//	}
	int p;
	for (p = 0; p < 4; p++) {
		fprintf(f, "SCALARS concentration_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (k = 0; k < I->gl_nz; k++) {
			for (i = 0; i < I->gl_nx; i++) {
				for (j = 0; j < I->gl_ny; j++) {
					if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
						fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, 4, i, j, k)]);
			}
		}
	}
//	fprintf(f, "COLOR_SCALARS saturation 3\n");
//	for (k = 0; k < I->gl_nz; k++) {
//		for (i = 0; i < I->gl_nx; i++) {
//			for (j = 0; j < I->gl_ny; j++) {
//				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
//					for (p = 5; p < 8; p++)
//						fprintf(f, "%20.10f\t", I->gl_B[GL_A_IND(I, p, i, j, k)]);
//					fprintf(f, "\n");
//				}
//			}
//		}
//	}
	for (p = 5; p < 8; p++) {
		fprintf(f, "SCALARS saturation_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (k = 0; k < I->gl_nz; k++) {
			for (i = 0; i < I->gl_nx; i++) {
				for (j = 0; j < I->gl_ny; j++) {
					if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
						fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS temperature_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, 8, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS temperature_environment double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, 9, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS processor int 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%d\n", I->ind_proc[i * I->gl_ny + j]);
			}
		}
	}
/*
//there is no function for avarage velocity
	fprintf(f, "VECTORS avarage_velocity_water double\n");
	for (k = 0; k < I->gl_nz; k++) {
		for (i = 0; i < I->gl_nx; i++) {
			for (j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.10f\t", avarage_velocity(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_oil double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.10f\t", avarage_velocity(I, 1, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_gas double\n");
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++)
						fprintf(f, "%20.10f\t", avarage_velocity(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
*/
	fclose(f);
	return 0;
}

int print_vtk(in *I, int n)
{
#if DEBUG
	printf("Process %d: printing result.\n", I->my_rank);
#endif
#if AVALANCHE
	if (print_vtk_avalanche(I, n)) return 1;
#endif
#if TERMOGAS
	if (I->nproc > 1) {
		//if (reconstruct_src(I)) return 1;
		if ((I->my_rank == 0) && (print_vtk_termogas_parallel(I, n))) return 1;
	} else {
		if (print_vtk_termogas(I, n)) return 1;
	}
#endif
	return 0;
}

int print_parameter_in_subdomains(in *I, int p, int n)
{
	int nn = n;
	int i = 0, j, k, a;
	while (nn > 0) {
		nn /= 10;
		i++;
	}
	//char file_name[8 + i];
	char file_name[50 + i];
	sprintf(file_name, "result/processor%d_parameter%d_time_step%d.log", I->my_rank, p, n);
	FILE *f = fopen(file_name, "w");
	for (k = 0; k < I->nz; k++) {
		fprintf(f, "k = %d\n", k);
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				fprintf(f, "%10lf\t", I->B_prev[B_IND(I, p, i, j, k)]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}
