#include "init_data.h"
#include "utils.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if AVALANCHE
int SET_initial_CONDITION_velocity_fixed_value(in *I)
{
	printf("Set the initial condition for velocity with fixed value for all calculation domain\n");
	int i, j, k, p;
	double fixed_value = 0;
	for (p = 0; p < 3; p++) {
		for (k = 0; k < I->nz; k++) {
			for (i = 0; i < I->nx; i++) {
				for (j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1) {
						I->B_prev[B_IND(I, p, i, j, k)] = fixed_value;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_phase_fraction_fixed_value(in *I)
{
	printf("Set the initial condition for velocity with fixed value for all calculation domain\n");
	int i, j, k;
	double fixed_value = I->density_snow;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					I->B_prev[B_IND(I, 3, i, j, k)] = fixed_value;
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_velocity_mass_is_mooving(in *I)
{
	printf("Set the initial condition for velocity with fixed value for all calculation domain\n");
	int i, j, k;
	double fixed_value[3];
	fixed_value[0] = 0.5;
	fixed_value[1] = 0;
	fixed_value[2] = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					if (phase_fraction(I, i, j, k) != 0) {
						I->B_prev[B_IND(I, 0, i, j, k)] = fixed_value[0];
						I->B_prev[B_IND(I, 1, i, j, k)] = fixed_value[1];
						I->B_prev[B_IND(I, 2, i, j, k)] = fixed_value[2];
					} else {
						I->B_prev[B_IND(I, 0, i, j, k)] = 0;
						I->B_prev[B_IND(I, 1, i, j, k)] = 0;
						I->B_prev[B_IND(I, 2, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_phase_fraction_ascii_map(in *I)
{
	printf("Set the initial condition for phase fraction with values token from the ASCII map\n");
	int i, j, k;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_boundary_cells[i * I->ny + j] != -1) {
					if (I->snow_region[(i / I->kx) * I->ncols + j / I->ky] == 1) {
						if ((k + 1) * (I->cellsize / (double) I->kz) <= I->depth) {
							I->B_prev[B_IND(I, 3, i, j, k)] = 1;
							I->mass_quantity += I->B_prev[B_IND(I, 3, i, j, k)];
						} else if (((k + 1) * (I->cellsize / (double) I->kz) > I->depth) && (k * (I->cellsize / (double) I->kz)) < I->depth) {
							I->B_prev[B_IND(I, 3, i, j, k)] = (I->depth - k * (I->cellsize / (double) I->kz)) / (I->cellsize / (double) I->kz);
							I->mass_quantity += I->B_prev[B_IND(I, 3, i, j, k)];
						} else {
							I->B_prev[B_IND(I, 3, i, j, k)] = 0;
						}
					} else {
						I->B_prev[B_IND(I, 3, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_pressure_fixed_value_with_hydrostatic_pressure(in *I)
{
	printf("Set the initial condition for pressure with fixed value atmosphere pressure in respect hydrostatic pressire\n");
	int i, j, k;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_boundary_cells[i * I->ny + j] != -1) {
					if (I->snow_region[(i / I->kx) * I->ncols + j / I->ky] == 1) {
						if ((k + 1) * (I->cellsize / (double) I->kz) <= I->depth) {
							I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere + I->density_snow * I->g[2] * (I->depth - k * (I->cellsize / (double) I->kz));
							//I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere;
						} else if (((k + 1) * (I->cellsize / (double) I->kz) > I->depth) && (k * (I->cellsize / (double) I->kz)) < I->depth) {
							I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere + I->density_snow * I->g[2] * (I->depth - k * (I->cellsize / (double) I->kz));
							//I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere;
						} else {
							I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere;
						}
					} else {
						I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere;
					}
				}
			}
		}
	}
	return 0;
}

int SET_initial_CONDITION_pressure_fixed_value(in *I)
{
	printf("Set the initial condition for pressure with fixed value atmosphere pressure in respect hydrostatic pressire\n");
	int i, j, k;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_boundary_cells[i * I->ny + j] != -1) {
					I->B_prev[B_IND(I, 4, i, j, k)] = I->pressure_atmosphere;
				}
			}
		}
	}
	return 0;
}
#endif

#if TERMOGAS
int SET_initial_CONDITION_termogas_fixed_value(in *I)
{
#if DEBUG
	printf("Set the initial condition for all parameters in termogas case with fixed value for all calculation domain in process %d\n", I->my_rank);
#endif
	for (int k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (int i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (int j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (internal_cell(I, i, j, k)) {
					if (well(I, i, j, k))
						I->porousness[POR_IND(I, i, j, k)] = I->porousness_value;
					else
						I->porousness[POR_IND(I, i, j, k)] = ((double) rand() / (RAND_MAX+1.0)) * 0.1 + I->porousness_value - 0.05;
					I->B_prev[B_IND(I, 0, i, j, k)] = 0;
					I->B_prev[B_IND(I, 1, i, j, k)] = 0;
					I->B_prev[B_IND(I, 2, i, j, k)] = 0;
					I->B_prev[B_IND(I, 3, i, j, k)] = 1;
					I->B_prev[B_IND(I, 4, i, j, k)] = I->initial_pressure;
					I->B_prev[B_IND(I, 5, i, j, k)] = 0;//I->residual_saturation[0];//2 * I->residual_saturation[0];// + I->epsilon;
					I->B_prev[B_IND(I, 6, i, j, k)] = 1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2];
					I->B_prev[B_IND(I, 7, i, j, k)] = 0;//I->residual_saturation[2];// + I->epsilon;
					I->B_prev[B_IND(I, 8, i, j, k)] = 0;//I->initial_temperature;
					I->B_prev[B_IND(I, 9, i, j, k)] = 0;//I->initial_temperature;
					if ((0) && (i + I->ind_start_region_proc[0] == I->gl_nx / 2) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2)) {
						I->B_prev[B_IND(I, 0, i, j, k)] = 3.55 * (1 - 2 * I->epsilon) / 4.55;
						for (int l = 0; l < 10; l++)
							for (int m = 0; m < 10; m++)
								I->B_prev[B_IND(I, 0, i - l, j - m, k)] = 1;//3.55 * (1 - 2 * I->epsilon) / 4.55;
						I->B_prev[B_IND(I, 1, i, j, k)] = (1 - 2 * I->epsilon) / 4.55;
						I->B_prev[B_IND(I, 2, i, j, k)] = I->epsilon;
						I->B_prev[B_IND(I, 3, i, j, k)] = I->epsilon;
						I->B_prev[B_IND(I, 4, i, j, k)] = I->injection_well_pressure;
						I->B_prev[B_IND(I, 5, i, j, k)] = 0.5 - (I->residual_saturation[1] + I->epsilon) / 2;
						I->B_prev[B_IND(I, 6, i, j, k)] = I->residual_saturation[1] + I->epsilon;
						I->B_prev[B_IND(I, 7, i, j, k)] = 0.5 - (I->residual_saturation[1] + I->epsilon) / 2;
						I->B_prev[B_IND(I, 8, i, j, k)] = I->injection_well_temperature - I->initial_temperature;
					}
				}
			}
		}
	}
	for (int k = -I->stencil_size; k < I->nz + I->stencil_size; k++)
		for (int i = -I->stencil_size; i < I->nx + I->stencil_size; i++)
			for (int j = -I->stencil_size; j < I->ny + I->stencil_size; j++)
				if (boundary_cell(I, i, j, k))
						I->porousness[POR_IND(I, i, j, k)] = I->porousness[POR_INT_NEI_IND(I, i, j, I->nz - 1)];
	return 0;
}
#endif
