#include "init_data.h"
#include "utils.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

int A_IND(in *I, int p, int i, int j, int k)
{
	return (I->num_parameters * (k * I->n_cells_multipl + I->ind_cell_multipl[i * I->ny + j]) + p);
}

int GL_A_IND(in *I, int p, int i, int j, int k)
{
#if DEBUG==4
	printf("I->gl_n_cells_multipl = %d\n", I->gl_n_cells_multipl);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_n_cells_multipl = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_n_cells_multipl);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_ind_cell_multipl[i * I->gl_ny + j] = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_ind_cell_multipl[i * I->gl_ny + j]);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_n_cells_multipl = %d, I->gl_ind_cell_multipl[i * I->gl_ny + j] = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_n_cells_multipl, I->gl_ind_cell_multipl[i * I->gl_ny + j]);
#endif
	return (I->num_parameters * (k * I->gl_n_cells_multipl + I->gl_ind_cell_multipl[i * I->gl_ny + j]) + p);
}

int B_IND(in *I, int p, int i, int j, int k)
{
	return (I->num_parameters * ((k + I->stencil_size) * I->n_boundary_cells + I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size]) + p);
}

int AREA_IND(in *I, int i, int j, int k, int s)
{
	return (6 * I->ind_cell_multipl[i * I->ny + j] + s);
}

int NORMAL_IND(in *I, int p, int i, int j, int k, int s)
{
	return (18 * I->ind_cell_multipl[i * I->ny + j] + s * 3 + p);
}

int VOLUME_IND(in *I, int i, int j, int k)
{
	return (I->ind_cell_multipl[i * I->ny + j]);
}

int check_for_corrupt_cell(in *I, int i, int j, int k)
{
	if ((i >= I->nx) || (i < 0) ||
		(j >= I->ny) || (j < 0) ||
		(k >= I->nz) || (k < 0) ||
		(I->ind_cell_multipl[i * I->ny + j] == -1)) {
			printf("!!!!!!Operating with corrupt cell!!!!!!\n");
			return 1;
	}
	return 0;
}

int well(in *I, int i, int j, int k)
{
	if ((injection_well(I, i, j, k)) || (production_well(I, i, j, k)))
		return 1;
	else
		return 0;
}

int production_well(in *I, int i, int j, int k)
{
		return 0;
	if ((i + I->ind_start_region_proc[0] == 0) && (j + I->ind_start_region_proc[1] == 0))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == 0) && (j + I->ind_start_region_proc[1] == I->gl_ny - 1))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx - 1) && (j + I->ind_start_region_proc[1] == 0))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx - 1) && (j + I->ind_start_region_proc[1] == I->gl_ny - 1))
		return 1;
	else
		return 0;
}

int injection_well(in *I, int i, int j, int k)
{
#if DEBUG==4
	printf("Process %d, PID %d: I->ind_start_region_proc[0] = %d\n", I->my_rank, getpid(), I->ind_start_region_proc[0]);
	printf("Process %d, PID %d: I->ind_start_region_proc[1] = %d\n", I->my_rank, getpid(), I->ind_start_region_proc[1]);
	printf("Process %d, PID %d: I->gl_nx = %d, I->nx = %d\n", I->my_rank, getpid(), I->gl_nx, I->nx);
	printf("Process %d, PID %d: I->gl_ny = %d, I->ny = %d\n", I->my_rank, getpid(), I->gl_ny, I->ny);
#endif
	if ((i + I->ind_start_region_proc[0] == I->gl_nx / 2) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2))
		return 1;
	else
		return 0;
}

int write_to_A_csr(in *I, int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value)
{
	if ((i >= 0) && (i < I->nx) && (j >= 0) && (j < I->ny) && (k >= 0) && (k < I->nz) && (I->ind_cell_multipl[i * I->ny + j] != -1) && (A_IND_S_SWITCH(I, i, j, k, s))) {
		int l, m, column, row;
		if (!((i_eq >= 0) && (i_eq < I->nx) && (j_eq >= 0) && (j_eq < I->ny) && (k_eq >= 0) && (k_eq < I->nz) && (I->ind_cell_multipl[i_eq * I->ny + j_eq] != -1)))
			printf("something very strange in ./src/lib/utils.c:87\n");
		if (s == -1) {
			column = A_IND(I, p, i, j, k);
		} else {
			column = A_IND_S(I, p, i, j, k, s);
		}
		row = A_IND(I, p_eq, i_eq, j_eq, k_eq);
		if (I->flag_first_time_step) {
			for (l = I->Aiptr_csr[row]; l <= I->A_ind_current; l++) {
				if (I->Ajptr_csr[l] == -1) {
					I->Aelem_csr[l] = value;
					I->Ajptr_csr[l] = column;
					break;
				} else if (I->Ajptr_csr[l] == column) {
					I->Aelem_csr[l] += value;
					break;
				} else if (I->Ajptr_csr[l] > column) {
					I->A_ind_current++;
					if (I->A_ind_current == I->non_zero_elem) {
						I->non_zero_elem++;
						if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					for (m = I->A_ind_current; m > l; m--) {
						I->Aelem_csr[m] = I->Aelem_csr[m - 1];
						I->Ajptr_csr[m] = I->Ajptr_csr[m - 1];
					}
					I->Aelem_csr[l] = value;
					I->Ajptr_csr[l] = column;
					break;
				} else if (l == I->A_ind_current) {
					I->A_ind_current++;
					if (I->A_ind_current == I->non_zero_elem) {
						I->non_zero_elem++;
						if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					I->Aelem_csr[I->A_ind_current] = value;
					I->Ajptr_csr[I->A_ind_current] = column;
					break;
				}
			}
		} else {
			for (l = I->Aiptr_csr[row]; l < I->Aiptr_csr[row + 1]; l++) {
				if (I->Ajptr_csr[l] == column) {
					I->Aelem_csr[l] += value;
					break;
				}
			}
		}
	} else {
		I->B[A_IND(I, p_eq, i_eq, j_eq, k_eq)] -= value * I->B_prev[B_IND(I, p, i, j, k)];
	}
	return 0;
}

int A_IND_S_SWITCH(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			if ((i + 1 >= I->nx) || (I->ind_cell_multipl[(i + 1) * I->ny + j] == -1))
				return 0;
			else
				return 1;
		case 1:
			if ((i - 1 < 0) || (I->ind_cell_multipl[(i - 1) * I->ny + j] == -1))
				return 0;
			else
				return 1;
		case 2:
			if ((j + 1 >= I->ny) || (I->ind_cell_multipl[i * I->ny + j + 1] == -1))
				return 0;
			else
				return 1;
		case 3:
			if ((j - 1 < 0) || (I->ind_cell_multipl[i * I->ny + j - 1] == -1))
				return 0;
			else
				return 1;
		case 4:
			if (k + 1 >= I->nz)
				return 0;
			else
				return 1;
		case 5:
			if (k - 1 < 0)
				return 0;
			else
				return 1;
	}
	return 1;
}

int A_IND_S(in *I, int p, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return A_IND(I, p, i + 1, j, k);
		case 1:
			return A_IND(I, p, i - 1, j, k);
		case 2:
			return A_IND(I, p, i, j + 1, k);
		case 3:
			return A_IND(I, p, i, j - 1, k);
		case 4:
			return A_IND(I, p, i, j, k + 1);
		case 5:
			return A_IND(I, p, i, j, k - 1);
	}
}

int cell_of_computation_domain(in *I, int i, int j, int k)
{
	if ((i >= - I->stencil_size) && (i < I->nx + I->stencil_size) &&
		(j >= - I->stencil_size) && (j < I->ny + I->stencil_size) &&
		(k >= - I->stencil_size) && (k < I->nz + I->stencil_size))
			return 1;
	else
			return 0;
}

int internal_cell(in *I, int i, int j, int k)
{
	if ((i >= 0) && (i < I->nx) && (j >= 0) && (j < I->ny) && (I->ind_cell_multipl[i * I->ny + j] != -1) && (k >= 0) && (k < I->nz))
		return 1;
	else
		return 0;

}

int boundary_cell(in *I, int i, int j, int k)
{
	if (cell_of_computation_domain(I, i, j, k) &&
		(I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
		(!(internal_cell(I, i, j, k))))
			return 1;
	else
			return 0;

}

int count_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 1, j, k)) x++;
	if (boundary_cell(I, i + 1, j, k)) x++;
	if (boundary_cell(I, i, j - 1, k)) x++;
	if (boundary_cell(I, i, j + 1, k)) x++;
	return x;
}

int count_second_order_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 2, j, k)) x++;
	if (boundary_cell(I, i + 2, j, k)) x++;
	if (boundary_cell(I, i, j - 2, k)) x++;
	if (boundary_cell(I, i, j + 2, k)) x++;
	return x;
}

int count_other_corner_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 1, j - 1, k)) x++;
	if (boundary_cell(I, i - 1, j + 1, k)) x++;
	if (boundary_cell(I, i + 1, j - 1, k)) x++;
	if (boundary_cell(I, i + 1, j + 1, k)) x++;
	if (boundary_cell(I, i - 1, j - 2, k)) x++;
	if (boundary_cell(I, i - 1, j + 2, k)) x++;
	if (boundary_cell(I, i + 1, j - 2, k)) x++;
	if (boundary_cell(I, i + 1, j + 2, k)) x++;
	if (boundary_cell(I, i - 2, j - 1, k)) x++;
	if (boundary_cell(I, i - 2, j + 1, k)) x++;
	if (boundary_cell(I, i + 2, j - 1, k)) x++;
	if (boundary_cell(I, i + 2, j + 1, k)) x++;
	if (boundary_cell(I, i - 2, j - 2, k)) x++;
	if (boundary_cell(I, i - 2, j + 2, k)) x++;
	if (boundary_cell(I, i + 2, j - 2, k)) x++;
	if (boundary_cell(I, i + 2, j + 2, k)) x++;
	return x;
}

#if AVALANCHE
double density(in *I, int i, int j, int k)
{
	//return phase_fraction(I, i, j, k) * I->density_snow + (1 - phase_fraction(I, i, j, k)) * I->density_air;
	return phase_fraction(I, i, j, k);
}

double velocity(in *I, int p, int i, int j, int k)
{
	return I->B_prev[B_IND(I, p, i, j, k)];
}

double phase_fraction(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 3, i, j, k)];
}

double pressure(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 4, i, j, k)];
}

double velocity_on_face(in *I, int p, int i, int j, int k, int s)
{
	int a, b;
	switch (s) {
		case 0:
			return (velocity(I, p, i, j, k) + velocity(I, p, i + 1, j, k)) / 2;
		case 1:
			return (velocity(I, p, i, j, k) + velocity(I, p, i - 1, j, k)) / 2;
		case 2:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j + 1, k)) / 2;
		case 3:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j - 1, k)) / 2;
		case 4:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j, k + 1)) / 2;
		case 5:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j, k - 1)) / 2;
	}
}

double density_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (density(I, i, j, k) + density(I, i + 1, j, k)) / 2;
		case 1:
			return (density(I, i, j, k) + density(I, i - 1, j, k)) / 2;
		case 2:
			return (density(I, i, j, k) + density(I, i, j + 1, k)) / 2;
		case 3:
			return (density(I, i, j, k) + density(I, i, j - 1, k)) / 2;
		case 4:
			return (density(I, i, j, k) + density(I, i, j, k + 1)) / 2;
		case 5:
			return (density(I, i, j, k) + density(I, i, j, k - 1)) / 2;
	}
}

double strain_rate_on_face(in *I, int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity_on_face(I, 0, i + 1, j, k, s) - velocity_on_face(I, 0, i - 1, j, k, s)) / (2 * I->dx[0]);
		case 1:
			return 0.5 * ((velocity_on_face(I, 0, i, j + 1, k, s) - velocity_on_face(I, 0, i, j - 1, k, s)) / (2 * I->dx[1]) +
					(velocity_on_face(I, 1, i + 1, j, k, s) - velocity_on_face(I, 1, i - 1, j, k, s)) / (2 * I->dx[0]));
		case 2:
			return 0.5 * ((velocity_on_face(I, 0, i, j, k + 1, s) - velocity_on_face(I, 0, i, j, k - 1, s)) / (2 * I->dx[2]) +
					(velocity_on_face(I, 2, i + 1, j, k, s) - velocity_on_face(I, 2, i - 1, j, k, s)) / (2 * I->dx[0]));
		case 3:
			return (velocity_on_face(I, 1, i, j + 1, k, s) - velocity_on_face(I, 1, i, j - 1, k, s)) / (2 * I->dx[1]);
		case 5:
			return 0.5 * ((velocity_on_face(I, 1, i, j, k + 1, s) - velocity_on_face(I, 1, i, j, k - 1, s)) / (2 * I->dx[2]) +
					(velocity_on_face(I, 2, i, j + 1, k, s) - velocity_on_face(I, 2, i, j - 1, k, s)) / (2 * I->dx[1]));
		case 8:
			return (velocity_on_face(I, 2, i, j, k + 1, s) - velocity_on_face(I, 2, i, j, k - 1, s)) / (2 * I->dx[2]);
	}
}

double shear_rate_on_face(in *I, int i, int j, int k, int s)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate_on_face(I, i, j, k, s, m, n) * strain_rate_on_face(I, i, j, k, s, m, n);
		}
	}
	return x;
}

double strain_rate(in *I, int i, int j, int k, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity(I, 0, i + 1, j, k) - velocity(I, 0, i - 1, j, k)) / (2 * I->dx[0]);
		case 1:
			return 0.5 * ((velocity(I, 0, i, j + 1, k) - velocity(I, 0, i, j - 1, k)) / (2 * I->dx[1]) +
				(velocity(I, 1, i + 1, j, k) - velocity(I, 1, i - 1, j, k)) / (2 * I->dx[0]));
		case 2:
			return 0.5 * ((velocity(I, 0, i, j, k + 1) - velocity(I, 0, i, j, k - 1)) / (2 * I->dx[2]) +
				(velocity(I, 2, i + 1, j, k) - velocity(I, 2, i - 1, j, k)) / (2 * I->dx[0]));
		case 3:
			return (velocity(I, 1, i, j + 1, k) - velocity(I, 1, i, j - 1, k)) / (2 * I->dx[1]);
		case 5:
			return 0.5 * ((velocity(I, 1, i, j, k + 1) - velocity(I, 1, i, j, k - 1)) / (2 * I->dx[2]) +
				(velocity(I, 2, i, j + 1, k) - velocity(I, 2, i, j - 1, k)) / (2 * I->dx[1]));
		case 8:
			return (velocity(I, 2, i, j, k + 1) - velocity(I, 2, i, j, k - 1)) / (2 * I->dx[2]);
	}
}

double shear_rate(in *I, int i, int j, int k)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate(I, i, j, k, m, n) * strain_rate(I, i, j, k, m, n);
		}
	}
	return x;
}

double phase_fraction_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i + 1, j, k)) / 2;
		case 1:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i - 1, j, k)) / 2;
		case 2:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j + 1, k)) / 2;
		case 3:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j - 1, k)) / 2;
		case 4:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j, k + 1)) / 2;
		case 5:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j, k - 1)) / 2;
	}
}

double effective_viscosity_on_face(in *I, int i, int j, int k, int s)
{
	double phase_fraction_on_face_tmp = phase_fraction_on_face(I, i, j, k, s);
	double shear_rate_on_face_tmp = shear_rate_on_face(I, i, j, k, s);
	double x = I->k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= I->shear_rate_0) {
		x += I->limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (I->k_viscosity_snow * pow(shear_rate_on_face_tmp, I->flow_index - 1) + I->yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

double pressure_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (pressure(I, i, j, k) + pressure(I, i + 1, j, k)) / 2;
		case 1:
			return (pressure(I, i, j, k) + pressure(I, i - 1, j, k)) / 2;
		case 2:
			return (pressure(I, i, j, k) + pressure(I, i, j + 1, k)) / 2;
		case 3:
			return (pressure(I, i, j, k) + pressure(I, i, j - 1, k)) / 2;
		case 4:
			return (pressure(I, i, j, k) + pressure(I, i, j, k + 1)) / 2;
		case 5:
			return (pressure(I, i, j, k) + pressure(I, i, j, k - 1)) / 2;
	}
}

int check_conservation_of_mass(in *I)
{
	int i, j, k;
	double mass = 0;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
					((!((i >= 0) && (i < I->nx) &&
					   (j >= 0) && (j < I->ny) && 
					   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
					 (k < 0) || (k >= I->nz))) {
							I->mass_quantity -= phase_fraction(I, i, j, k);
				}
			}
		}
	}
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					mass += phase_fraction(I, i, j, k);
				}
			}
		}
	}
	if (mass == I->mass_quantity) {
		return 0;
	} else {
		return 1;
	}
}

int barotropy_pressure(in *I)
{
	int i, j, k;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					if (check_for_corrupt_cell(I, i, j, k)) return 1;
					I->B[A_IND(I, 4, i, j, k)] = I->pressure_atmosphere * (
						1 +
						0.001 * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) / I->density_snow +
						0 * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) / (I->density_snow * I->density_snow));
				}
			}
		}
	}
	return 0;
}

int barotropy_density(in *I)
{
	int i, j, k;
	double c1 = 1, c2 = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					if (check_for_corrupt_cell(I, i, j, k)) return 1;
					I->B[A_IND(I, 3, i, j, k)] = (I->density_snow / I->pressure_atmosphere) * (
						I->pressure_atmosphere +
						c1 * (I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere) +
						c2 * 0.5 * (I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere) *
						(I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere));
				}
			}
		}
	}
	return 0;
}
#endif

#if TERMOGAS
double saturation(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2))) {
		printf("Error saturation index\n");
		return 0;
	}
	return I->B_prev[B_IND(I, p + 5, i, j, k)];
}

double concentration(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error concentration index\n");
		return 0;
	}
	return I->B_prev[B_IND(I, p, i, j, k)];
}

double pressure(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 4, i, j, k)];
}

double temperature_flow(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 8, i, j, k)];
}

double temperature_environment(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 9, i, j, k)];
}

double density_t(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1))
		return ((I->density_0[p] + pow(I->density_coef_a[p], -2) * (pressure(I, i, j, k) - I->pressure_0)) / (1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0)));
	else if (p == 2)
		return (pressure(I, i, j, k) * (concentration(I, 0, i, j, k) * I->molar_weight[0] +
										concentration(I, 1, i, j, k) * I->molar_weight[1] +
										concentration(I, 2, i, j, k) * I->molar_weight[2] +
										concentration(I, 3, i, j, k) * I->molar_weight[3]) / (I->R * temperature_flow(I, i, j, k)));
	else {
		printf("Error density index\n");
		return 0;
	}
}

double two_phase_relative_permeability(in *I, int p, int pr, int i, int j, int k)
{
	double S;
	if ((p == 0) && (pr == 1)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]) / (1 - I->residual_saturation[0] - I->residual_saturation[2]);
		return (pow(S, 4));
	}
	else if ((p == 1) && (pr == 0)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]) / (1 - I->residual_saturation[0] - I->residual_saturation[2]);
		return ((1 - S) * (1 - S) * (1 - S * S));
	}
	else if ((p == 0) && (pr == 2)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]) / (1 - I->residual_saturation[0] - I->residual_saturation[1]);
		return (pow(S, 4));
	}
	else if ((p == 2) && (pr == 0)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]) / (1 - I->residual_saturation[0] - I->residual_saturation[1]);
		return ((1 - S) * (1 - S) * (1 - S * S));
	}
	else if ((p == 1) && (pr == 2)) {
		S = (saturation(I, 1, i, j, k) - I->residual_saturation[1]) / (1 - I->residual_saturation[1] - I->residual_saturation[0]);
		return (pow(S, 4));
	}
	else if ((p == 2) && (pr == 1)) {
		S = (saturation(I, 1, i, j, k) - I->residual_saturation[1]) / (1 - I->residual_saturation[1] - I->residual_saturation[0]);
		return ((1 - S) * (1 - S) * (1 - S * S));
	} else {
		printf("Error two phase relative permeability index\n");
		return 0;
	}
}

double relative_permeability(in *I, int p, int i, int j, int k)
{
	if (p == 0)
		return ((saturation(I, 1, i, j, k) - I->residual_saturation[1]) * two_phase_relative_permeability(I, 0, 1, i, j, k) +
				(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * two_phase_relative_permeability(I, 0, 2, i, j, k)) /
				(saturation(I, 1, i, j, k) - I->residual_saturation[1] + saturation(I, 2, i, j, k) - I->residual_saturation[2]);
	else if (p == 1)
		return ((saturation(I, 0, i, j, k) - I->residual_saturation[0]) * two_phase_relative_permeability(I, 1, 0, i, j, k) +
				(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * two_phase_relative_permeability(I, 1, 2, i, j, k)) /
				(saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 2, i, j, k) - I->residual_saturation[2]);
	else if (p == 2)
		return ((saturation(I, 0, i, j, k) - I->residual_saturation[0]) * two_phase_relative_permeability(I, 2, 0, i, j, k) +
				(saturation(I, 1, i, j, k) - I->residual_saturation[1]) * two_phase_relative_permeability(I, 2, 1, i, j, k)) /
				(saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 1, i, j, k) - I->residual_saturation[1]);
	else {
		printf("Error relative permeability index\n");
		return 0;
	}
}

double viscosity_gas(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error viscosity gas index\n");
		return 0;
	}
#if DEBUG==3
	printf("Process %d: I->viscosity_coef_A_gas[%d] = %lf\n", I->my_rank, p, I->viscosity_coef_A_gas[p]);
	printf("Process %d: I->temperature_0_gas[%d] = %lf\n", I->my_rank, p, I->temperature_0_gas[p]);
	printf("Process %d: I->viscosity_coef_C_gas[%d] = %lf\n", I->my_rank, p, I->viscosity_coef_C_gas[p]);
	printf("Process %d: temperature_flow(I, %d, %d, %d) = %lf\n", I->my_rank, i, j, k, temperature_flow(I, i, j, k));
#endif
	//printf("Process %d: \n", I->my_rank);
	double tmp =  I->viscosity_coef_A_gas[p] *
		(I->temperature_0_gas[p] + I->viscosity_coef_C_gas[p]) *
		pow(temperature_flow(I, i, j, k) / I->temperature_0_gas[p], 3 / 2) /
		(temperature_flow(I, i, j, k) + I->viscosity_coef_C_gas[p]);
	if (tmp < 0)
		printf("Process %d, PID %d: p= %d, i = %d, j = %d, k = %d, viscosity_gas = %lf\n", I->my_rank, getpid(), p, i, j, k, tmp);
	return tmp;
}

double molar_fraction(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error molar fraction index\n");
		return 0;
	}
	double x = 0;
	int l;
	for (l = 0; l < 4; l++)
		x += concentration(I, l, i, j, k) / I->molar_weight[l];
	return concentration(I, p, i, j, k) / (I->molar_weight[p] * x);
}

double viscosity_water(in *I, int i, int j, int k)
{
	double tmp;
	//tmp = 0.001 / (0.14 + (temperature_flow(I, i, j, k) - 273.15) / 30 + 0.000009 * pow(temperature_flow(I, i, j, k) - 273.15, 2));
	//tmp = (I->viscosity_coef_A[0] / (1 / density_t(I, 0, i, j, k) - I->viscosity_coef_B[0]));
	tmp = 0.001 / (0.14 + (temperature_flow(I, i, j, k) - 273.15) / 30.0 + 0.000009 * (temperature_flow(I, i, j, k) - 273.15) * (temperature_flow(I, i, j, k) - 273.15));
	if (tmp < 0)
		printf("Process %d, PID %d: i = %d, j = %d, k = %d, viscosity_water = %lf, temperature_flow = %lf\n", I->my_rank, getpid(), i, j, k, tmp, temperature_flow(I, i, j, k));
	return tmp;
}

double viscosity_oil(in *I, int i, int j, int k)
{
	double tmp;
	//tmp = density_t(I, 1, i, j, k) * (pow(0.00000236, pow(I->initial_temperature / temperature_flow(I, i, j, k), 4)) - 0.0000006);
	//tmp = (I->viscosity_coef_A[1] / (1 / density_t(I, 1, i, j, k) - I->viscosity_coef_B[1]));
	tmp = density_t(I, 1, i, j, k) * 0.000001 * (pow((1000.0 + 0.6), pow((30.0 / (temperature_flow(I, i, j, k) - 273.15)), 4.0)) - 0.6);
	if (tmp < 0)
		printf("Process %d, PID %d: i = %d, j = %d, k = %d, viscosity_oil = %lf, temperature_flow = %lf\n", I->my_rank, getpid(), i, j, k, tmp, temperature_flow(I, i, j, k));
	return tmp;
}

double viscosity(in *I, int p, int i, int j, int k)
{
	if (p == 0) {
		//return (I->viscosity_coef_A[p] / (1 / density_t(I, p, i, j, k) - I->viscosity_coef_B[p]));
		return viscosity_water(I, i, j, k);
	} else if (p == 1) {
		//return (I->viscosity_coef_A[p] / (1 / density_t(I, p, i, j, k) - I->viscosity_coef_B[p]));
		return viscosity_oil(I, i, j, k);
	} else if (p == 2) {
		double x = 1;
		int l;
		for (l = 0; l < 4; l++) {
#if DEBUG==3
			printf("Process %d: l = %d, i = %d, j = %d, k = %d\n", I->my_rank, l, i, j, k);
			printf("Process %d: viscosity_gas(I, l, i, j, k) = %lf\n", I->my_rank, viscosity_gas(I, l, i, j, k));
			printf("Process %d: molar_fraction(I, l, i, j, k) = %lf\n", I->my_rank, molar_fraction(I, l, i, j, k));
#endif
			x *= pow(viscosity_gas(I, l, i, j, k), molar_fraction(I, l, i, j, k));
		}
		return x;
	} else {
		printf("Error viscosity index\n");
		return 0;
	}
}

double Darsi_M_coef_phases(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2))) {
		printf("Error Darsi M coefficient phases index\n");
		return 0;
	}
	return density_t(I, p, i, j, k) * I->permeability * relative_permeability(I, p, i, j, k) / viscosity(I, p, i, j, k);
}

double Darsi_M_coef(in *I, int i, int j, int k)
{
	double x = 0;
	int l;
	for (l = 0; l < 3; l++)
		x += Darsi_M_coef_phases(I, l, i, j, k);
	//return 1;
	return x;
}

double capillary_pressure_derivative_by_saturation(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 2))) {
		printf("Error capillary pressure index\n");
		return 0;
	}
	return - I->capillary_pressure_at_maximum_saturation[p] *
		pow(1 - I->residual_saturation_two_phase[p], I->capillary_pressure_coef) *
		I->capillary_pressure_coef *
		pow(saturation(I, p, i, j, k), - I->capillary_pressure_coef - 1);
}

double grad_pressure(in *I, int pr, int i, int j, int k)
{
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	return (pressure(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - pressure(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
}

double grad_saturation(in *I, int p, int pr, int i, int j, int k)
{
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	return (saturation(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - saturation(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
}

double avarage_velocity(in *I, int p, int pr, int i, int j, int k) //p - oil, water, gas; pr - x1, x2, x3
{
	if (!((pr == 0) || (pr == 1) || (pr == 2))) {
		printf("Error avarage velocity index\n");
		return 0;
	}
	if (pr != 2) {
		return 1;
	} else {
		return 0;
	}
	printf("WTF?!\n");
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	if (p == 0) {
#if DEBUG
		if (viscosity(I, p, i, j, k) < 0.0000001) {
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: viscosity = %lf\n", I->my_rank, getpid(), p, i, j, k, viscosity(I, p, i, j, k));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: density_t(I, p, i, j, k) = %lf\n", I->my_rank, getpid(), p, i, j, k, density_t(I, p, i, j, k));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: numerator of density = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, (I->density_0[p] + pow(I->density_coef_a[p], -2) * (pressure(I, i, j, k) - I->pressure_0)));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: denominator of density = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, (1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0)));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: temperature_flow = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, temperature_flow(I, i, j, k));
		}
#endif
		return - I->permeability * relative_permeability(I, p, i, j, k) *(
			grad_pressure(I, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k) - 1) *
			capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *
			grad_saturation(I, 0, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k)) *
			capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *
			grad_saturation(I, 2, pr, i, j, k)
			) / viscosity(I, p, i, j, k);
	} else if (p == 1) {
#if DEBUG
		if (viscosity(I, p, i, j, k) < 0.0000001) {
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: viscosity = %lf\n", I->my_rank, getpid(), p, i, j, k, viscosity(I, p, i, j, k));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: density_t(I, p, i, j, k) = %lf\n", I->my_rank, getpid(), p, i, j, k, density_t(I, p, i, j, k));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: numerator of density = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, (I->density_0[p] + pow(I->density_coef_a[p], -2) * (pressure(I, i, j, k) - I->pressure_0)));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: denominator of density = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, (1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0)));
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: temperature_flow = %lf\n",\
				I->my_rank, getpid(), p, i, j, k, temperature_flow(I, i, j, k));
		}
#endif
		return - I->permeability * relative_permeability(I, p, i, j, k) * (
			grad_pressure(I, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k)) *
			capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *
			grad_saturation(I, 0, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k)) *
			capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *
			grad_saturation(I, 2, pr, i, j, k)
			) / viscosity(I, p, i, j, k);
	} else if (p == 2) {
#if DEBUG
		if (viscosity(I, p, i, j, k) < 0.0000001)
			printf("Process %d, PID %d, p = %d, i = %d, j = %d, k = %d: viscosity = %lf\n", I->my_rank, getpid(), p, i, j, k, viscosity(I, p, i, j, k));
#endif
		return - I->permeability * relative_permeability(I, p, i, j, k) * (
			grad_pressure(I, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k)) *
			capillary_pressure_derivative_by_saturation(I, 0, i, j, k) *
			grad_saturation(I, 0, pr, i, j, k) +
			(Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k) - 1) *
			capillary_pressure_derivative_by_saturation(I, 2, i, j, k) *
			grad_saturation(I, 2, pr, i, j, k)
			) / viscosity(I, p, i, j, k);
	} else {
		printf("Error avarage velocity index\n");
		return 0;
	}
}

double rate_of_reaction_coef(in *I, int i, int j, int k)
{
	if (temperature_flow(I, i, j, k) < I->threshold_temperature)
		return 0;
	else
		return 0;
//		return I->frequency_factor * exp(- I->activation_temperature / temperature_flow(I, i, j, k));
//	if (injection_well(I, i, j, k))
//		return 0;
//	else
}

double rate_of_reaction(in *I, int i, int j, int k)
{
//	double tmp = rate_of_reaction_coef(I, i, j, k) * concentration(I, 1, i, j, k) * pow(relative_permeability(I, 1, i, j, k) * relative_permeability(I, 2, i, j, k), 0.25);
//	double tmp = rate_of_reaction_coef(I, i, j, k) * concentration(I, 1, i, j, k) * saturation(I, 2, i, j, k);
//	double tmp = rate_of_reaction_coef(I, i, j, k) * (saturation(I, 1, i, j, k) * density_t(I, 1, i, j, k) / I->molar_weight[4]) *
//		pow((concentration(I, 1, i, j, k) * saturation(I, 2, i, j, k) *
//			(pressure(I, i, j, k) * I->molar_weight[1] / (I->R * temperature_flow(I, i, j, k))) / I->molar_weight[1]), 12.5);
//	double tmp = rate_of_reaction_coef(I, i, j, k) * (saturation(I, 1, i, j, k) * density_t(I, 1, i, j, k) / I->molar_weight[4]) *
//		(concentration(I, 1, i, j, k) * saturation(I, 2, i, j, k) *
//		(pressure(I, i, j, k) * I->molar_weight[1] / (I->R * temperature_flow(I, i, j, k))) / I->molar_weight[1]);
	double tmp = I->porousness * rate_of_reaction_coef(I, i, j, k) * (saturation(I, 1, i, j, k) - I->residual_saturation[1]) *
		(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * concentration(I, 1, i, j, k) *
		pow(pressure(I, i, j, k) / I->pressure_activ, I->stoichiometric_coef_activ);
	return tmp;
	//return 0;
}

double rate_of_reaction_derivative_by_temperature(in *I, int i, int j, int k)
{
	double tmp = rate_of_reaction(I, i, j, k) * I->activation_temperature / (temperature_flow(I, i, j, k) * temperature_flow(I, i, j, k));
	return tmp;
	//return 0;
}

double mass_inflow_rate_func(in *I, int p, int i, int j, int k)
{
	if (p == 0) // N2
		return 0;
	else if (p == 1) // O2
		return rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 2) // CO2
		return rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 3) // H2O(g)
		return rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 5) // water
		return rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[3];
	else if (p == 6) // oil
		return rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[4];
	else if (p == 7) // gas
		return mass_inflow_rate_func(I, 0, i, j, k) + mass_inflow_rate_func(I, 1, i, j, k) +
			mass_inflow_rate_func(I, 2, i, j, k) + mass_inflow_rate_func(I, 3, i, j, k);
	else {
		printf("Error mass inflow rate index\n");
		return 0;
	}
}

double mass_inflow_rate_func_derivative_by_temperature(in *I, int p, int i, int j, int k)
{
	if (p == 0) // N2
		return 0;
	else if (p == 1) // O2
		return rate_of_reaction_derivative_by_temperature(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 2) // CO2
		return rate_of_reaction_derivative_by_temperature(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 3) // H2O
		return rate_of_reaction_derivative_by_temperature(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 5) // water
		return rate_of_reaction_derivative_by_temperature(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[3];
	else if (p == 6) // oil
		return rate_of_reaction_derivative_by_temperature(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[4];
	else if (p == 7) // gas
		return mass_inflow_rate_func_derivative_by_temperature(I, 0, i, j, k) + mass_inflow_rate_func_derivative_by_temperature(I, 1, i, j, k) +
			mass_inflow_rate_func_derivative_by_temperature(I, 2, i, j, k) + mass_inflow_rate_func_derivative_by_temperature(I, 3, i, j, k);
	else {
		printf("Error mass inflow rate index\n");
		return 0;
	}
}

double chemical_reaction_heat_flow(in *I, int i, int j, int k)
{
	int p;
	double x = 0;
	for (p = 0; p < 4; p++) {
		x -= mass_inflow_rate_func(I, p, i, j, k) * I->initial_enthalpy[p + 2];
	}
	for (p = 0; p < 2; p++) {
		x -= mass_inflow_rate_func(I, p + 5, i, j, k) * I->initial_enthalpy[p];
	}
	return x;
}

double chemical_reaction_heat_flow_derivative_by_temperature(in *I, int i, int j, int k)
{
	return I->heat_liberation * mass_inflow_rate_func_derivative_by_temperature(I, 6, i, j, k);
}

double density_derivative_by_pressure(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1))
		return pow(I->density_coef_a[p], -2) / (1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0));
	else if (p == 2)
		return (concentration(I, 0, i, j, k) * I->molar_weight[0] +
				concentration(I, 1, i, j, k) * I->molar_weight[1] +
				concentration(I, 2, i, j, k) * I->molar_weight[2] +
				concentration(I, 3, i, j, k) * I->molar_weight[3]) / (I->R * temperature_flow(I, i, j, k));
	else {
		printf("Error density derivative by pressure index\n");
		return 0;
	}
}

double Darsi_A_coef(in *I, int i, int j, int k)
{
	int l;
	double x = 0;
	for (l = 0; l < 3; l++) {
		x += I->porousness * saturation(I, l, i, j, k) * density_derivative_by_pressure(I, l, i, j, k);
	}
	return x;
}

double internal_energy(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1)) // water, oil
		return I->specific_heat[p] * (temperature_flow(I, i, j, k) - I->tempetarure_for_calculation_internal_energy) + I->initial_enthalpy[p];
	else if (p == 2) // gas
		return internal_energy(I, 3, i, j, k) + internal_energy(I, 4, i, j, k) + internal_energy(I, 5, i, j, k) + internal_energy(I, 6, i, j, k);
	else if ((p == 3) || (p == 4) || (p == 5) || (p == 6)) // N2, O2, CO2, H2O
		return (I->specific_heat[p - 1] * (temperature_flow(I, i, j, k) - I->tempetarure_for_calculation_internal_energy) + I->initial_enthalpy[p - 1]) * concentration(I, p - 3, i, j, k);
	else if (p == 7) // environment
		return I->specific_heat[p - 1] * (temperature_environment(I, i, j, k) - I->tempetarure_for_calculation_internal_energy) + I->initial_enthalpy[p - 1];
	else {
		printf("Error internal energy index\n");
		return 0;
	}
}

double enthalpy_flow(in *I, int i, int j, int k)
{
	return internal_energy(I, 0, i, j, k) - internal_energy(I, 1, i, j, k) - internal_energy(I, 2, i, j, k) + 2 * internal_energy(I, 5, i, j, k) + 2 * internal_energy(I, 6, i, j, k);
}

int check_sum(in *I)
{
	int i, j, k, p;
	double C = 0;
	for (p = 0; p < 3; p++)
		C += I->residual_saturation[p] + I->epsilon;
	double saturation_sum, concentration_sum;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++) {
						if (saturation(I, p, i, j, k) < I->residual_saturation[p] + I->epsilon) {
							//printf("Too small saturation of %d phase at [%d, %d, %d]: %.20f\n", p, i, j, k, saturation(I, p, i, j, k));
							I->B_prev[B_IND(I, p + 5, i, j, k)] = I->residual_saturation[p] + I->epsilon;
						}
					}
					saturation_sum = 0;
					for (p = 0; p < 3; p++)
						saturation_sum += saturation(I, p, i, j, k);
					if (saturation_sum != 1) {
						//printf("Saturatins sum is not equal 1 at [%d, %d, %d]: %.20f\n", i, j, k, saturation_sum);
						I->B_prev[B_IND(I, 5, i, j, k)] = I->residual_saturation[0] + I->epsilon + (saturation(I, 0, i, j, k) - I->residual_saturation[0] - I->epsilon) * (1 - C) / (saturation_sum - C);
						I->B_prev[B_IND(I, 6, i, j, k)] = I->residual_saturation[1] + I->epsilon + (saturation(I, 1, i, j, k) - I->residual_saturation[1] - I->epsilon) * (1 - C) / (saturation_sum - C);
						I->B_prev[B_IND(I, 7, i, j, k)] = I->residual_saturation[2] + I->epsilon + (saturation(I, 2, i, j, k) - I->residual_saturation[2] - I->epsilon) * (1 - C) / (saturation_sum - C);
					}
					for (p = 0; p < 4; p++) {
						if (concentration(I, p, i, j, k) < I->epsilon) {
							//printf("Too small concentration of %d component at [%d, %d, %d]: %.20f\n", p, i, j, k, concentration(I, p, i, j, k));
							I->B_prev[B_IND(I, p, i, j, k)] = I->epsilon;
						}
					}
					concentration_sum = 0;
					for (p = 0; p < 4; p++)
						concentration_sum += concentration(I, p, i, j, k);
					if (concentration(I, 0, i, j, k) + concentration(I, 1, i, j, k) + concentration(I, 2, i, j, k) + concentration(I, 3, i, j, k) != 1) {
						//printf("Concentration sum is not equal 1 at [%d, %d, %d]: %.20f\n", i, j, k, concentration(I, 0, i, j, k) + concentration(I, 1, i, j, k) + concentration(I, 2, i, j, k) + concentration(I, 3, i, j, k));
						I->B_prev[B_IND(I, 0, i, j, k)] = I->epsilon + (concentration(I, 0, i, j, k) - I->epsilon) * (1 - 4 * I->epsilon) / concentration_sum;
						I->B_prev[B_IND(I, 1, i, j, k)] = I->epsilon + (concentration(I, 1, i, j, k) - I->epsilon) * (1 - 4 * I->epsilon) / concentration_sum;
						I->B_prev[B_IND(I, 2, i, j, k)] = I->epsilon + (concentration(I, 2, i, j, k) - I->epsilon) * (1 - 4 * I->epsilon) / concentration_sum;
						I->B_prev[B_IND(I, 3, i, j, k)] = I->epsilon + (concentration(I, 3, i, j, k) - I->epsilon) * (1 - 4 * I->epsilon) / concentration_sum;
					}
				}
			}
		}
	}
	return 0;
}

int print_oil_production(in *I)
{
	FILE *f;
	int well_coordinates[3];
	well_coordinates[0] = well_coordinates[1] = well_coordinates[2] = 0;
	double velocity_of_oil_production;
	/*well_coordinates[0] = I->nx / 2;
	well_coordinates[1] = I->ny / 2;
	well_coordinates[2] = 0;*/
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_oil_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf", I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness * density_t(I, 1, 1, 1, 1));
	} else {
		if ((f = fopen("result/velocity_of_oil_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	velocity_of_oil_production = saturation(I, 1, well_coordinates[0], well_coordinates[1], well_coordinates[2]) * I->porousness * I->dx[0] * I->dx[1] * I->dx[2] -
		I->B[A_IND(I, 6, well_coordinates[0], well_coordinates[1], well_coordinates[2])] * I->porousness * I->dx[0] * I->dx[1] * I->dx[2] -
		avarage_velocity(I, 1, 0, well_coordinates[0] + 1, well_coordinates[1], well_coordinates[2]) * I->dx[1] * I->dx[2] * I->dt -
		avarage_velocity(I, 1, 1, well_coordinates[0], well_coordinates[1] + 1, well_coordinates[2]) * I->dx[0] * I->dx[2] * I->dt;
	fprintf(f, "%lf\t%lf\n", (I->time_step + 1) * I->dt, velocity_of_oil_production);
	fclose(f);

	I->volume_producted_oil_m += velocity_of_oil_production;
	if (I->time_step == 0) {
		if ((f = fopen("result/oil_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf\n%lf\t%lf\n", I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness * density_t(I, 1, 1, 1, 1), 0.0, 0.0);
	} else {
		if ((f = fopen("result/oil_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", (I->time_step + 1) * I->dt, I->volume_producted_oil_m);
	fclose(f);

	velocity_of_oil_production *= density_t(I, 1, well_coordinates[0], well_coordinates[1], well_coordinates[2]);
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_oil_production_kg.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf", I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness * density_t(I, 1, 1, 1, 1));
	} else {
		if ((f = fopen("result/velocity_of_oil_production_kg.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", (I->time_step + 1) * I->dt, velocity_of_oil_production);
	fclose(f);

	I->volume_producted_oil_kg += velocity_of_oil_production;
	if (I->time_step == 0) {
		if ((f = fopen("result/oil_production_kg.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf\n%lf\t%lf\n", I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness * density_t(I, 1, 1, 1, 1), 0.0, 0.0);
	} else {
		if ((f = fopen("result/oil_production_kg.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", (I->time_step + 1) * I->dt, I->volume_producted_oil_kg);
	fclose(f);
	return 0;
}

double avarage_velocity_global(in *I)
{
	int i, j, k, p, pr;
	double x = 0, y = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++) {
					//	for (pr = 0; pr < 3; pr++) {
					//		x += abs(avarage_velocity(I, p, pr, i, j, k));
					//	}
						x = sqrt(pow(avarage_velocity(I, p, 0, i, j, k), 2) + pow(avarage_velocity(I, p, 1, i, j, k), 2) + pow(avarage_velocity(I, p, 2, i, j, k), 2));
						if (y < x)
							y = x;
					}
				}
			}
		}
	}
//	x /= I->n_cells_multipl * 9;
	return y;
}

#endif
