#include "init_data.h"
#include "utils.h"
#include "vtk_map_functions.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int print_vtk_avalanche(in *I)
{
	char file_name[14 + 20];
	sprintf(file_name, "result/map%d.vtk", (int) I->time);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "slope\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
//	printf("n_points_multipl = %d\n", I->n_points_multipl);
//	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	fprintf(f, "POINTS %d double\n", I->n_points_multipl * ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	for (int k = 0; k <= (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 0; i < (I->nrows - 1) * I->kx + 1; i ++) {
			for (int j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
				if (I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * I->cellsize / (double) I->kx, j * I->cellsize / (double) I->ky,
						I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] + k * I->cellsize / (double) I->kz);
			}
		}
	}
	fprintf(f, "CELLS %d %d\n", I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize),
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) * 9);
	for (int k = 0; k < (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
			for (int j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
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
	for (int i = 0; i < I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize); i++) {
		fprintf(f, "%d\n", 12);
	}
	fprintf(f, "CELL_DATA %d\n", I->n_cells_multipl * I->nz);
	fprintf(f, "VECTORS velocity double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\t%20.10f\t%20.10f\n", I->B_prev[B_IND(I, 0, i, j, k)], I->B_prev[B_IND(I, 1, i, j, k)], I->B_prev[B_IND(I, 2, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\n", I->B_prev[B_IND(I, 4, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS snow_volume_fraction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.10f\n", I->B_prev[B_IND(I, 3, i, j, k)]);
			}
		}
	}
	fclose(f);
	sprintf(file_name, "tmp/B_prev%10.10lf.txt", I->time);
	f = fopen(file_name, "w");
	for (int p = 0; p < I->num_parameters; p++) {
		fprintf(f, "parameter %d\n", p);
		for (int k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			fprintf(f, "k = %d\n", k);
			for (int i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (int j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
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

int print_vtk_termogas(in *I)
{
	char file_name[14 + 20];
	int file_num;
	if (I->units_for_write_interval == 0) {
		if (I->time_step % I->write_interval)
			return 0;
		file_num = I->time_step / I->write_interval;
	} else {
		if (((int) I->time / I->units_for_write_interval) / I->write_interval == I->written_step)
			return 0;
		file_num = I->written_step = ((int) I->time / I->units_for_write_interval) / I->write_interval;
	}
	sprintf(file_name, "result/map%d.vtk", file_num);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
//	fprintf(f, "#time = %lf\n", I->time);
	fprintf(f, "oil displacement, time = %lf\n", I->time);
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
//	printf("n_points_multipl = %d\n", I->n_points_multipl);
//	printf("points in z direction %d\n", ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	fprintf(f, "POINTS %d double\n", I->n_points_multipl * ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1));
	for (int k = 0; k <= (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 0; i < (I->nrows - 1) * I->kx + 1; i ++) {
			for (int j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
				if (I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * I->cellsize / (double) I->kx, j * I->cellsize / (double) I->ky,
						I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] + k * I->cellsize / (double) I->kz);
			}
		}
	}
	fprintf(f, "CELLS %d %d\n", I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize),
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) * 9);
	for (int k = 0; k < (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
			for (int j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
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
	for (int i = 0; i < I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize); i++) {
		fprintf(f, "%d\n", 12);
	}
	fprintf(f, "CELL_DATA %d\n", I->n_cells_multipl * I->nz);
/*
	tmp = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					tmp += I->B_prev[B_IND(I, 0, i, j, k)];
			}
		}
	}
	printf("Summ = %lf\n", tmp);
*/
	for (int p = 0; p < 4; p++) {
		fprintf(f, "SCALARS concentration_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (int k = 0; k < I->nz; k++) {
			for (int i = 0; i < I->nx; i++) {
				for (int j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						fprintf(f, "%20.20f\n", concentration(I, p, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", pressure(I, i, j, k));
			}
		}
	}
	for (int p = 0; p < 3; p++) {
		fprintf(f, "SCALARS saturation_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (int k = 0; k < I->nz; k++) {
			for (int i = 0; i < I->nx; i++) {
				for (int j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						fprintf(f, "%20.20f\n", saturation(I, p, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS temperature_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", temperature_flow(I, i, j, k));
			}
		}
	}
	fprintf(f, "SCALARS temperature_environment double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1)
					fprintf(f, "%20.20f\n", temperature_environment(I, i, j, k));
			}
		}
	}
	for (int p = 0; p < 4; p++) {
		fprintf(f, "SCALARS concentration_%d_saturation_gas double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (int k = 0; k < I->nz; k++) {
			for (int i = 0; i < I->nx; i++) {
				for (int j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						fprintf(f, "%20.20f\n", concentration(I, p, i, j, k) * saturation(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS rate_of_reaction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", rate_of_reaction(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas_N2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity_gas_component(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas_O2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity_gas_component(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas_CO2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity_gas_component(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS viscosity_gas_H2O(g) double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", viscosity_gas_component(I, 3, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS saturation_water_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 5, i, j, k) /\
						(density_t(I, 0, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)]));
				}
			}
		}
	}
	fprintf(f, "SCALARS saturation_oil_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 6, i, j, k) /\
						(density_t(I, 1, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)]));
				}
			}
		}
	}
	fprintf(f, "SCALARS saturation_gas_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 7, i, j, k) /\
						(density_t(I, 2, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)]));
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_N2_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 0, i, j, k) /\
						(density_t(I, 2, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)] * saturation(I, 2, i, j, k)));
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_O2_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 1, i, j, k) /\
						(density_t(I, 2, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)] * saturation(I, 2, i, j, k)));
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_CO2_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 2, i, j, k) /\
						(density_t(I, 2, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)] * saturation(I, 2, i, j, k)));
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_H2O(g)_inflow_rate double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n",\
						mass_inflow_rate(I, 3, i, j, k) /\
						(density_t(I, 2, i, j, k) * I->dx[0] * I->dx[1] * I->dx[2] * I->porousness[POR_IND(I, i, j, k)] * saturation(I, 2, i, j, k)));
				}
			}
		}
	}
/*
	fprintf(f, "SCALARS molar_fraction_N2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", molar_fraction(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS molar_fraction_O2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", molar_fraction(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS molar_fraction_CO2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", molar_fraction(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS molar_fraction_H2O(g) double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", molar_fraction(I, 3, i, j, k));
				}
			}
		}
	}
	*/
	fprintf(f, "SCALARS density_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS density_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS density_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", density_t(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS rate_of_reaction_coef double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", rate_of_reaction_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS pressure_in_rate_of_reaction double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", pow(pressure(I, i, j, k) / I->pressure_activ, I->stoichiometric_coef_activ));
				}
			}
		}
	}
	fprintf(f, "SCALARS sat_res_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 0, i, j, k) - I->residual_saturation[0]);
				}
			}
		}
	}
	fprintf(f, "SCALARS sat_res_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 1, i, j, k) - I->residual_saturation[1]);
				}
			}
		}
	}
	fprintf(f, "SCALARS sat_res_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 2, i, j, k) - I->residual_saturation[2]);
				}
			}
		}
	}
	/*
	fprintf(f, "SCALARS press_oxigen double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", pow(pressure(I, i, j, k) / I->pressure_activ, I->stoichiometric_coef_activ));
				}
			}
		}
	}
	fprintf(f, "SCALARS rate_coef double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					double tmp = I->porousness * (saturation(I, 1, i, j, k) - I->residual_saturation[1]) *
						(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * concentration(I, 1, i, j, k) *
						pow(pressure(I, i, j, k) / I->pressure_activ, I->stoichiometric_coef_activ);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					double tmp = (concentration(I, 1, i, j, k) * saturation(I, 2, i, j, k) *
						(pressure(I, i, j, k) * I->molar_weight[1] / (I->R * temperature_flow(I, i, j, k))) / I->molar_weight[1]);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS concentration_oxigen double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					double tmp = (saturation(I, 1, i, j, k) * density_t(I, 1, i, j, k) / I->molar_weight[4]);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	*/
	fprintf(f, "SCALARS chemical_reaction_heat_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", chemical_reaction_heat_flow(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS chemical_reaction_heat_flow_K_per_sec double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					double tmp = 0;
					for (int p = 0; p < 2; p++)
						tmp += I->porousness[POR_IND(I, i, j, k)] * density_t(I, p, i, j, k) * saturation(I, p, i, j, k) * I->specific_heat[p];
					for (int p = 0; p < 4; p++)
						tmp += I->porousness[POR_IND(I, i, j, k)] * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[p + 2] * concentration(I, p, i, j, k);
					fprintf(f, "%20.20f\n", chemical_reaction_heat_flow(I, i, j, k) / tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 0, i, j, k) * mass_inflow_rate(I, 5, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 1, i, j, k) * mass_inflow_rate(I, 6, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", saturation(I, 2, i, j, k) * mass_inflow_rate(I, 7, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_N2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", concentration(I, 0, i, j, k) * mass_inflow_rate(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_O2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", concentration(I, 1, i, j, k) * mass_inflow_rate(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_CO2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", concentration(I, 2, i, j, k) * mass_inflow_rate(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_H2O double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", concentration(I, 3, i, j, k) * mass_inflow_rate(I, 3, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS mass_inflow_rate_summ double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					double tmp =	saturation(I, 0, i, j, k) * mass_inflow_rate(I, 5, i, j, k) +
							saturation(I, 1, i, j, k) * mass_inflow_rate(I, 6, i, j, k) +
							saturation(I, 2, i, j, k) * mass_inflow_rate(I, 7, i, j, k);
					fprintf(f, "%20.20f\n", tmp);
				}
			}
		}
	}
	fprintf(f, "SCALARS porousness double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", I->porousness[POR_IND(I, i, j, k)]);
				}
			}
		}
	}
	/*
	fprintf(f, "SCALARS avarage_velocity_coef_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS avarage_velocity_coef_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef_water_rel double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef_gas_rel double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS Darsi_M_coef double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", Darsi_M_coef(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS capillary_pressure_derivative_by_saturation_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", capillary_pressure_derivative_by_saturation(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS capillary_pressure_derivative_by_saturation_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", capillary_pressure_derivative_by_saturation(I, 2, i, j, k));
				}
			}
		}
	}
	*/
	fprintf(f, "SCALARS relative_permeability_water double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", relative_permeability(I, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS relative_permeability_oil double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", relative_permeability(I, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS relative_permeability_gas double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", relative_permeability(I, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_0 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 0, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_1 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 1, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_2 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 0, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_3 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 2, 0, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_4 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 1, 2, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS two_phase_relative_permeability_5 double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", two_phase_relative_permeability(I, 2, 1, i, j, k));
				}
			}
		}
	}
	fprintf(f, "SCALARS thermal_conductivity double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					fprintf(f, "%20.20f\n", thermal_conductivity(I, i, j, k));
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_water double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 0, p, i, j, k) / (I->porousness[POR_IND(I, i, j, k)] * saturation(I, 0, i, j, k)));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS filtration_rate_oil double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 1, p, i, j, k) / (I->porousness[POR_IND(I, i, j, k)] * saturation(I, 1, i, j, k)));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_oil double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 1, p, i, j, k) / (I->porousness[POR_IND(I, i, j, k)] * saturation(I, 1, i, j, k)));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_gas double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 2, p, i, j, k) / (I->porousness[POR_IND(I, i, j, k)] * saturation(I, 2, i, j, k)));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_water_divide_by_saturation double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 0, p, i, j, k) / (saturation(I, 0, i, j, k) * I->porousness[POR_IND(I, i, j, k)]));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_oil_divide_by_saturation double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 1, p, i, j, k) / (saturation(I, 1, i, j, k) * I->porousness[POR_IND(I, i, j, k)]));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS avarage_velocity_gas_divide_by_saturation double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 2, p, i, j, k) / (saturation(I, 2, i, j, k) * I->porousness[POR_IND(I, i, j, k)]));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS grad_pressure double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", grad_pressure(I, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	/*
	fprintf(f, "VECTORS vel_water_of_sat_water double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k) - 1) *\
							capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *\
							grad_saturation(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS vel_water_of_sat_gas double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k)) *\
							capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *\
							grad_saturation(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS vel_oil_of_sat_water double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k)) *\
							capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *\
							grad_saturation(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS vel_oil_of_sat_gas double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k)) *\
							capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *\
							grad_saturation(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS vel_gas_of_sat_water double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k)) *\
							capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *\
							grad_saturation(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS vel_gas_of_sat_gas double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t",\
							(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k) - 1) *\
							capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *\
							grad_saturation(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS grad_saturation_water double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", grad_saturation(I, 0, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS grad_saturation_oil double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", grad_saturation(I, 1, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS grad_saturation_gas double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", grad_saturation(I, 2, p, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}
	fprintf(f, "VECTORS velosity_concentration double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int p = 0; p < 3; p++)
						fprintf(f, "%20.20f\t", avarage_velocity(I, 2, p, i, j, k) / (saturation(I, 2, i, j, k) * I->porousness));
					fprintf(f, "\n");
				}
			}
		}
	}
	*/
	fprintf(f, "VECTORS velosity_temperature_flow double\n");
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (int pr = 0; pr < 3; pr++)
						fprintf(f, "%20.20f\t", convection_speed(I, pr, i, j, k));
					fprintf(f, "\n");
				}
			}
		}
	}


	fclose(f);
	return 0;
}

int print_vtk_termogas_parallel(in *I)
{
	char file_name[14 + 20];
	sprintf(file_name, "result/map%d.vtk", I->time_step);
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
	for (int k = 0; k <= (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 0; i < (I->nrows - 1) * I->kx + 1; i ++) {
			for (int j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
				if (I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1)
					fprintf(f, "%f %f %f\n", i * I->cellsize / (double) I->kx, j * I->cellsize / (double) I->ky,
						I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] + k * I->cellsize / (double) I->kz);
			}
		}
	}
	fprintf(f, "CELLS %d %d\n", I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize),
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) * 9);
	for (int k = 0; k < (int) (I->hight / I->cellsize) * I->kz; k++) {
		for (int i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
			for (int j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
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
	for (int i = 0; i < I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize); i++) {
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
	for (int p = 0; p < 4; p++) {
		fprintf(f, "SCALARS concentration_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (int k = 0; k < I->gl_nz; k++) {
			for (int i = 0; i < I->gl_nx; i++) {
				for (int j = 0; j < I->gl_ny; j++) {
					if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
						fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS pressure double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->gl_nz; k++) {
		for (int i = 0; i < I->gl_nx; i++) {
			for (int j = 0; j < I->gl_ny; j++) {
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
	for (int p = 5; p < 8; p++) {
		fprintf(f, "SCALARS saturation_%d double 1\n", p);
		fprintf(f, "LOOKUP_TABLE default\n");
		for (int k = 0; k < I->gl_nz; k++) {
			for (int i = 0; i < I->gl_nx; i++) {
				for (int j = 0; j < I->gl_ny; j++) {
					if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
						fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, p, i, j, k)]);
				}
			}
		}
	}
	fprintf(f, "SCALARS temperature_flow double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->gl_nz; k++) {
		for (int i = 0; i < I->gl_nx; i++) {
			for (int j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, 8, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS temperature_environment double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->gl_nz; k++) {
		for (int i = 0; i < I->gl_nx; i++) {
			for (int j = 0; j < I->gl_ny; j++) {
				if (I->gl_ind_cell_multipl[i * I->gl_ny + j] != -1)
					fprintf(f, "%20.20f\n", I->gl_B[GL_A_IND(I, 9, i, j, k)]);
			}
		}
	}
	fprintf(f, "SCALARS processor int 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int k = 0; k < I->gl_nz; k++) {
		for (int i = 0; i < I->gl_nx; i++) {
			for (int j = 0; j < I->gl_ny; j++) {
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

int print_vtk(in *I)
{
#if DEBUG
	printf("Process %d: printing result.\n", I->my_rank);
#endif
#if AVALANCHE
	if (print_vtk_avalanche(I)) return 1;
#endif
#if TERMOGAS
	if (I->nproc > 1) {
		//if (reconstruct_src(I)) return 1;
		if ((I->my_rank == 0) && (print_vtk_termogas_parallel(I))) return 1;
	} else {
		if (print_vtk_termogas(I)) return 1;
	}
#endif
	return 0;
}

int print_parameter_in_subdomains(in *I, int p)
{
	char file_name[50 + 20];
	sprintf(file_name, "result/processor%d_parameter%d_time%10.10lf.log", I->my_rank, p, I->time);
	FILE *f = fopen(file_name, "w");
	for (int k = 0; k < I->nz; k++) {
		fprintf(f, "k = %d\n", k);
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				fprintf(f, "%10lf\t", I->B_prev[B_IND(I, p, i, j, k)]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}
