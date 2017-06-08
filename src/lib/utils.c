#include "init_data.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int A_IND(in I, int p, int i, int j, int k)
{
	return (I.num_parameters * (k * I.n_cells_multipl + I.ind_cell_multipl[i * I.ny + j]) + p);
}

int B_IND(in I, int p, int i, int j, int k)
{
	return (I.num_parameters * ((k + I.stencil_size) * I.n_boundary_cells + I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size]) + p);
}

int AREA_IND(in I, int i, int j, int k, int s)
{
	return (6 * I.ind_cell_multipl[i * I.ny + j] + s);
}

int NORMAL_IND(in I, int p, int i, int j, int k, int s)
{
	return (18 * I.ind_cell_multipl[i * I.ny + j] + s * 3 + p);
}

int VOLUME_IND(in I, int i, int j, int k)
{
	return (I.ind_cell_multipl[i * I.ny + j]);
}

double density(in I, int i, int j, int k)
{
	return phase_fraction(I, i, j, k) * I.density_snow + (1 - phase_fraction(I, i, j, k)) * I.density_air;
}

double velocity(in I, int p, int i, int j, int k)
{
	return I.B_prev[B_IND(I, p, i, j, k)];
}

double phase_fraction(in I, int i, int j, int k)
{
	return I.B_prev[B_IND(I, 3, i, j, k)];
}

double pressure(in I, int i, int j, int k)
{
	return I.B_prev[B_IND(I, 4, i, j, k)];
}

int check_for_corrupt_cell(in I, int i, int j, int k)
{
	if ((i >= I.nx) || (i < 0) ||
		(j >= I.ny) || (j < 0) ||
		(k >= I.nz) || (k < 0) ||
		(I.ind_cell_multipl[i * I.ny + j] == -1)) {
			printf("!!!!!!Operating with corrupt cell!!!!!!\n");
			return 1;
	}
	return 0;
}

int write_to_A_csr(in I, int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value)
{
	if ((i >= 0) && (i < I.nx) && (j >= 0) && (j < I.ny) && (k >= 0) && (k < I.nz) && (I.ind_cell_multipl[i * I.ny + j] != -1) && (A_IND_S_SWITCH(I, i, j, k, s))) {
		int l, m, column, row;
		if (!((i_eq >= 0) && (i_eq < I.nx) && (j_eq >= 0) && (j_eq < I.ny) && (k_eq >= 0) && (k_eq < I.nz) && (I.ind_cell_multipl[i_eq * I.ny + j_eq] != -1)))
			printf("something wery strange\n");
		if (s == -1) {
			column = A_IND(I, p, i, j, k);
		} else {
			column = A_IND_S(I, p, i, j, k, s);
		}
		row = A_IND(I, p_eq, i_eq, j_eq, k_eq);
		if (I.flag_first_time_step) {
			for (l = I.Aiptr_csr[row]; l <= I.A_ind_current; l++) {
				if (I.Ajptr_csr[l] == -1) {
					I.Aelem_csr[l] = value;
					I.Ajptr_csr[l] = column;
					break;
				} else if (I.Ajptr_csr[l] == column) {
					I.Aelem_csr[l] += value;
					break;
				} else if (I.Ajptr_csr[l] > column) {
					I.A_ind_current++;
					if (I.A_ind_current == I.non_zero_elem) {
						I.non_zero_elem++;
						if ((I.Aelem_csr = (double *) realloc(I.Aelem_csr, I.non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I.Ajptr_csr = (int *) realloc(I.Ajptr_csr, I.non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					for (m = I.A_ind_current; m > l; m--) {
						I.Aelem_csr[m] = I.Aelem_csr[m - 1];
						I.Ajptr_csr[m] = I.Ajptr_csr[m - 1];
					}
					I.Aelem_csr[l] = value;
					I.Ajptr_csr[l] = column;
					break;
				} else if (l == I.A_ind_current) {
					I.A_ind_current++;
					if (I.A_ind_current == I.non_zero_elem) {
						I.non_zero_elem++;
						if ((I.Aelem_csr = (double *) realloc(I.Aelem_csr, I.non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I.Ajptr_csr = (int *) realloc(I.Ajptr_csr, I.non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					I.Aelem_csr[I.A_ind_current] = value;
					I.Ajptr_csr[I.A_ind_current] = column;
					break;
				}
			}
		} else {
			for (l = I.Aiptr_csr[row]; l < I.Aiptr_csr[row + 1]; l++) {
				if (I.Ajptr_csr[l] == column) {
					I.Aelem_csr[l] += value;
					break;
				}
			}
		}
	} else {
		I.B[A_IND(I, p_eq, i_eq, j_eq, k_eq)] -= value * I.B_prev[B_IND(I, p, i, j, k)];
	}
	return 0;
}

int A_IND_S_SWITCH(in I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			if ((i + 1 >= I.nx) || (I.ind_cell_multipl[(i + 1) * I.ny + j] == -1))
				return 0;
			else
				return 1;
		case 1:
			if ((i - 1 < 0) || (I.ind_cell_multipl[(i - 1) * I.ny + j] == -1))
				return 0;
			else
				return 1;
		case 2:
			if ((j + 1 >= I.ny) || (I.ind_cell_multipl[i * I.ny + j + 1] == -1))
				return 0;
			else
				return 1;
		case 3:
			if ((j - 1 < 0) || (I.ind_cell_multipl[i * I.ny + j - 1] == -1))
				return 0;
			else
				return 1;
		case 4:
			if (k + 1 >= I.nz)
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

int A_IND_S(in I, int p, int i, int j, int k, int s)
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

double velocity_on_face(in I, int p, int i, int j, int k, int s)
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

double density_on_face(in I, int i, int j, int k, int s)
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

double strain_rate_on_face(in I, int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity_on_face(I, 0, i + 1, j, k, s) - velocity_on_face(I, 0, i - 1, j, k, s)) / (2 * I.dx[0]);
		case 1:
			return 0.5 * ((velocity_on_face(I, 0, i, j + 1, k, s) - velocity_on_face(I, 0, i, j - 1, k, s)) / (2 * I.dx[1]) +
					(velocity_on_face(I, 1, i + 1, j, k, s) - velocity_on_face(I, 1, i - 1, j, k, s)) / (2 * I.dx[0]));
		case 2:
			return 0.5 * ((velocity_on_face(I, 0, i, j, k + 1, s) - velocity_on_face(I, 0, i, j, k - 1, s)) / (2 * I.dx[2]) +
					(velocity_on_face(I, 2, i + 1, j, k, s) - velocity_on_face(I, 2, i - 1, j, k, s)) / (2 * I.dx[0]));
		case 3:
			return (velocity_on_face(I, 1, i, j + 1, k, s) - velocity_on_face(I, 1, i, j - 1, k, s)) / (2 * I.dx[1]);
		case 5:
			return 0.5 * ((velocity_on_face(I, 1, i, j, k + 1, s) - velocity_on_face(I, 1, i, j, k - 1, s)) / (2 * I.dx[2]) +
					(velocity_on_face(I, 2, i, j + 1, k, s) - velocity_on_face(I, 2, i, j - 1, k, s)) / (2 * I.dx[1]));
		case 8:
			return (velocity_on_face(I, 2, i, j, k + 1, s) - velocity_on_face(I, 2, i, j, k - 1, s)) / (2 * I.dx[2]);
	}
}

double shear_rate_on_face(in I, int i, int j, int k, int s)
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

double phase_fraction_on_face(in I, int i, int j, int k, int s)
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

double effective_viscosity_on_face(in I, int i, int j, int k, int s)
{
	double phase_fraction_on_face_tmp = phase_fraction_on_face(I, i, j, k, s);
	double shear_rate_on_face_tmp = shear_rate_on_face(I, i, j, k, s);
	double x = I.k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= I.shear_rate_0) {
		x += I.limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (I.k_viscosity_snow * pow(shear_rate_on_face_tmp, I.flow_index - 1) + I.yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

double pressure_on_face(in I, int i, int j, int k, int s)
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

int check_conservation_of_mass(in I)
{
	int i, j, k;
	double mass = 0;
	for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
		for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
			for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
				if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
					((!((i >= 0) && (i < I.nx) &&
					   (j >= 0) && (j < I.ny) && 
					   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
					 (k < 0) || (k >= I.nz))) {
							I.mass_quantity -= phase_fraction(I, i, j, k);
				}
			}
		}
	}
	for (k = 0; k < I.nz; k++) {
		for (i = 0; i < I.nx; i++) {
			for (j = 0; j < I.ny; j++) {
				if (I.ind_cell_multipl[i * I.ny + j] != -1) {
					mass += phase_fraction(I, i, j, k);
				}
			}
		}
	}
	if (mass == I.mass_quantity) {
		return 0;
	} else {
		return 1;
	}
}

