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
#if AVALANCHE
	free(I->snow_region);
#endif
	free(I->Aelem_csr);
	free(I->Ajptr_csr);
	free(I->Aiptr_csr);
	free(I->B_prev);
	free(I->ind_boundary_cells);
	if (I->my_rank == 0) {
		free(I->ind_multipl);
		if (free_initial_arrays(I)) return 1;
	}
	if (I->nproc > 1) {
		if (free_parallel_arrays(I)) return 1;
	}
	free(I->ind_start_region_proc);
	free(I->array_of_parameters);
	return 0;
}

int free_parallel_arrays(in *I)
{
	free(I->gl_ind_cell_multipl);
	free(I->ind_proc);
	free(I->gl_B);
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

int set_array_of_parameters_termogas(in *I)
{
	if ((I->flag_first_time_step) && ((I->array_of_parameters = (double *) malloc(I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * I->dependent_variables * sizeof(double))) == NULL)) {
		printf("Memory error\n");
		return 1;
	}
	for (int k = - I->stencil_size + 1; k < I->gl_nz + I->stencil_size - 1; k++) {
		for (int i = - I->stencil_size + 1; i < I->gl_nx + I->stencil_size - 1; i++) {
			for (int j = - I-> stencil_size + 1; j < I->gl_ny + I->stencil_size - 1; j++) {
				I->array_of_parameters[PARAM_IND(I, 0, i, j, k)] = density_t_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 1, i, j, k)] = density_t_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 2, i, j, k)] = density_t_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 3, i, j, k)] = two_phase_relative_permeability_func(I, 0, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 4, i, j, k)] = two_phase_relative_permeability_func(I, 1, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 5, i, j, k)] = two_phase_relative_permeability_func(I, 0, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 6, i, j, k)] = two_phase_relative_permeability_func(I, 2, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 7, i, j, k)] = two_phase_relative_permeability_func(I, 1, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 8, i, j, k)] = two_phase_relative_permeability_func(I, 2, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 9, i, j, k)] = relative_permeability_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 10, i, j, k)] = relative_permeability_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 11, i, j, k)] = relative_permeability_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 12, i, j, k)] = viscosity_gas_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 13, i, j, k)] = viscosity_gas_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 14, i, j, k)] = viscosity_gas_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 15, i, j, k)] = viscosity_gas_func(I, 3, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 16, i, j, k)] = molar_fraction_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 17, i, j, k)] = molar_fraction_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 18, i, j, k)] = molar_fraction_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 19, i, j, k)] = molar_fraction_func(I, 3, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 20, i, j, k)] = viscosity_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 21, i, j, k)] = viscosity_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 22, i, j, k)] = viscosity_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 23, i, j, k)] = Darsi_M_coef_phases_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 24, i, j, k)] = Darsi_M_coef_phases_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 25, i, j, k)] = Darsi_M_coef_phases_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 26, i, j, k)] = Darsi_M_coef_func(I, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 27, i, j, k)] = capillary_pressure_derivative_by_saturation_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 28, i, j, k)] = capillary_pressure_derivative_by_saturation_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 29, i, j, k)] = grad_pressure_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 30, i, j, k)] = grad_pressure_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 31, i, j, k)] = grad_pressure_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 32, i, j, k)] = grad_saturation_func(I, 0, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 33, i, j, k)] = grad_saturation_func(I, 0, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 34, i, j, k)] = grad_saturation_func(I, 0, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 35, i, j, k)] = grad_saturation_func(I, 1, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 36, i, j, k)] = grad_saturation_func(I, 1, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 37, i, j, k)] = grad_saturation_func(I, 1, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 38, i, j, k)] = grad_saturation_func(I, 2, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 39, i, j, k)] = grad_saturation_func(I, 2, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 40, i, j, k)] = grad_saturation_func(I, 2, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 41, i, j, k)] = coef_grad_saturation_func(I, 0, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 42, i, j, k)] = coef_grad_saturation_func(I, 0, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 43, i, j, k)] = coef_grad_saturation_func(I, 1, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 44, i, j, k)] = coef_grad_saturation_func(I, 1, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 45, i, j, k)] = coef_grad_saturation_func(I, 2, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 46, i, j, k)] = coef_grad_saturation_func(I, 2, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 47, i, j, k)] = avarage_velocity_func(I, 0, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 48, i, j, k)] = avarage_velocity_func(I, 0, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 49, i, j, k)] = avarage_velocity_func(I, 0, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 50, i, j, k)] = avarage_velocity_func(I, 1, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 51, i, j, k)] = avarage_velocity_func(I, 1, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 52, i, j, k)] = avarage_velocity_func(I, 1, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 53, i, j, k)] = avarage_velocity_func(I, 2, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 54, i, j, k)] = avarage_velocity_func(I, 2, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 55, i, j, k)] = avarage_velocity_func(I, 2, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 56, i, j, k)] = rate_of_reaction_coef_func(I, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 57, i, j, k)] = rate_of_reaction_func(I, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 58, i, j, k)] = mass_inflow_rate_func(I, 0, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 59, i, j, k)] = mass_inflow_rate_func(I, 1, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 60, i, j, k)] = mass_inflow_rate_func(I, 2, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 61, i, j, k)] = mass_inflow_rate_func(I, 3, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 62, i, j, k)] = mass_inflow_rate_func(I, 5, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 63, i, j, k)] = mass_inflow_rate_func(I, 6, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 64, i, j, k)] = mass_inflow_rate_func(I, 7, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 65, i, j, k)] = chemical_reaction_heat_flow_func(I, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 66, i, j, k)] = Darsi_A_coef_func(I, i, j, k);
				I->array_of_parameters[PARAM_IND(I, 67, i, j, k)] = thermal_conductivity_func(I, i, j, k);
			}
		}
	}
	return 0;
}
