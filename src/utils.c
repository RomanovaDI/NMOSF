#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

int A_IND(int p, int i, int j, int k)
{
	return (num_parameters * (k * n_cells_multipl + ind_cell_multipl[i * ny + j]) + p);
}

int B_IND(int p, int i, int j, int k)
{
	return (num_parameters * ((k + stencil_size) * n_boundary_cells + ind_boundary_cells[(i + stencil_size) * (ny + 2 * stencil_size) + j + stencil_size]) + p);
}

int AREA_IND(int i, int j, int k, int s)
{
	return (6 * ind_cell_multipl[i * ny + j] + s);
}

int NORMAL_IND(int p, int i, int j, int k, int s)
{
	return (18 * ind_cell_multipl[i * ny + j] + s * 3 + p);
}

int VOLUME_IND(int i, int j, int k)
{
	return (ind_cell_multipl[i * ny + j]);
}

double density(int i, int j, int k)
{
	return phase_fraction(i, j, k) * density_snow + (1 - phase_fraction(i, j, k)) * density_air;
}

double velocity(int p, int i, int j, int k)
{
	return B_prev[B_IND(p, i, j, k)];
}

double phase_fraction(int i, int j, int k)
{
	return B_prev[B_IND(3, i, j, k)];
}

double pressure(int i, int j, int k)
{
	return B_prev[B_IND(4, i, j, k)];
}

int check_for_corrupt_cell(int i, int j, int k)
{
	if ((i >= nx) || (i < 0) ||
		(j >= ny) || (j < 0) ||
		(k >= nz) || (k < 0) ||
		(ind_cell_multipl[i * ny + j] == -1)) {
			printf("!!!!!!Operating with corrupt cell!!!!!!\n");
			return 1;
	}
	return 0;
}

int write_to_A_csr(int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value)
{
	if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz) && (ind_cell_multipl[i * ny + j] != -1) && (A_IND_S_SWITCH(i, j, k, s))) {
		int l, m, column, row;
		if (!((i_eq >= 0) && (i_eq < nx) && (j_eq >= 0) && (j_eq < ny) && (k_eq >= 0) && (k_eq < nz) && (ind_cell_multipl[i_eq * ny + j_eq] != -1)))
			printf("something wery strange\n");
		if (s == -1) {
			column = A_IND(p, i, j, k);
		} else {
			column = A_IND_S(p, i, j, k, s);
		}
		row = A_IND(p_eq, i_eq, j_eq, k_eq);
		if (flag_first_time_step) {
			for (l = Aiptr_csr[row]; l <= A_ind_current; l++) {
				if (Ajptr_csr[l] == -1) {
					Aelem_csr[l] = value;
					Ajptr_csr[l] = column;
					break;
				} else if (Ajptr_csr[l] == column) {
					Aelem_csr[l] += value;
					break;
				} else if (Ajptr_csr[l] > column) {
					A_ind_current++;
					if (A_ind_current == non_zero_elem) {
						non_zero_elem++;
						if ((Aelem_csr = (double *) realloc(Aelem_csr, non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((Ajptr_csr = (int *) realloc(Ajptr_csr, non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					for (m = A_ind_current; m > l; m--) {
						Aelem_csr[m] = Aelem_csr[m - 1];
						Ajptr_csr[m] = Ajptr_csr[m - 1];
					}
					Aelem_csr[l] = value;
					Ajptr_csr[l] = column;
					break;
				} else if (l == A_ind_current) {
					A_ind_current++;
					if (A_ind_current == non_zero_elem) {
						non_zero_elem++;
						if ((Aelem_csr = (double *) realloc(Aelem_csr, non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((Ajptr_csr = (int *) realloc(Ajptr_csr, non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					Aelem_csr[A_ind_current] = value;
					Ajptr_csr[A_ind_current] = column;
					break;
				}
			}
		} else {
			for (l = Aiptr_csr[row]; l < Aiptr_csr[row + 1]; l++) {
				if (Ajptr_csr[l] == column) {
					Aelem_csr[l] += value;
					break;
				}
			}
		}
	} else {
		B[A_IND(p_eq, i_eq, j_eq, k_eq)] -= value * B_prev[B_IND(p, i, j, k)];
	}
	return 0;
}

int A_IND_S_SWITCH(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			if ((i + 1 >= nx) || (ind_cell_multipl[(i + 1) * ny + j] == -1))
				return 0;
			else
				return 1;
		case 1:
			if ((i - 1 < 0) || (ind_cell_multipl[(i - 1) * ny + j] == -1))
				return 0;
			else
				return 1;
		case 2:
			if ((j + 1 >= ny) || (ind_cell_multipl[i * ny + j + 1] == -1))
				return 0;
			else
				return 1;
		case 3:
			if ((j - 1 < 0) || (ind_cell_multipl[i * ny + j - 1] == -1))
				return 0;
			else
				return 1;
		case 4:
			if (k + 1 >= nz)
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

int A_IND_S(int p, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return A_IND(p, i + 1, j, k);
		case 1:
			return A_IND(p, i - 1, j, k);
		case 2:
			return A_IND(p, i, j + 1, k);
		case 3:
			return A_IND(p, i, j - 1, k);
		case 4:
			return A_IND(p, i, j, k + 1);
		case 5:
			return A_IND(p, i, j, k - 1);
	}
}

double velocity_on_face(int p, int i, int j, int k, int s)
{
	int a, b;
	switch (s) {
		case 0:
			return (velocity(p, i, j, k) + velocity(p, i + 1, j, k)) / 2;
		case 1:
			return (velocity(p, i, j, k) + velocity(p, i - 1, j, k)) / 2;
		case 2:
			return (velocity(p, i, j, k) + velocity(p, i, j + 1, k)) / 2;
		case 3:
			return (velocity(p, i, j, k) + velocity(p, i, j - 1, k)) / 2;
		case 4:
			return (velocity(p, i, j, k) + velocity(p, i, j, k + 1)) / 2;
		case 5:
			return (velocity(p, i, j, k) + velocity(p, i, j, k - 1)) / 2;
	}
}

double density_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (density(i, j, k) + density(i + 1, j, k)) / 2;
		case 1:
			return (density(i, j, k) + density(i - 1, j, k)) / 2;
		case 2:
			return (density(i, j, k) + density(i, j + 1, k)) / 2;
		case 3:
			return (density(i, j, k) + density(i, j - 1, k)) / 2;
		case 4:
			return (density(i, j, k) + density(i, j, k + 1)) / 2;
		case 5:
			return (density(i, j, k) + density(i, j, k - 1)) / 2;
	}
}

double strain_rate_on_face(int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity_on_face(0, i + 1, j, k, s) - velocity_on_face(0, i - 1, j, k, s)) / (2 * dx[0]);
		case 1:
			return 0.5 * ((velocity_on_face(0, i, j + 1, k, s) - velocity_on_face(0, i, j - 1, k, s)) / (2 * dx[1]) +
					(velocity_on_face(1, i + 1, j, k, s) - velocity_on_face(1, i - 1, j, k, s)) / (2 * dx[0]));
		case 2:
			return 0.5 * ((velocity_on_face(0, i, j, k + 1, s) - velocity_on_face(0, i, j, k - 1, s)) / (2 * dx[2]) +
					(velocity_on_face(2, i + 1, j, k, s) - velocity_on_face(2, i - 1, j, k, s)) / (2 * dx[0]));
		case 3:
			return (velocity_on_face(1, i, j + 1, k, s) - velocity_on_face(1, i, j - 1, k, s)) / (2 * dx[1]);
		case 5:
			return 0.5 * ((velocity_on_face(1, i, j, k + 1, s) - velocity_on_face(1, i, j, k - 1, s)) / (2 * dx[2]) +
					(velocity_on_face(2, i, j + 1, k, s) - velocity_on_face(2, i, j - 1, k, s)) / (2 * dx[1]));
		case 8:
			return (velocity_on_face(2, i, j, k + 1, s) - velocity_on_face(2, i, j, k - 1, s)) / (2 * dx[2]);
	}
}

double shear_rate_on_face(int i, int j, int k, int s)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate_on_face(i, j, k, s, m, n) * strain_rate_on_face(i, j, k, s, m, n);
		}
	}
	return x;
}

double phase_fraction_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (phase_fraction(i, j, k) + phase_fraction(i + 1, j, k)) / 2;
		case 1:
			return (phase_fraction(i, j, k) + phase_fraction(i - 1, j, k)) / 2;
		case 2:
			return (phase_fraction(i, j, k) + phase_fraction(i, j + 1, k)) / 2;
		case 3:
			return (phase_fraction(i, j, k) + phase_fraction(i, j - 1, k)) / 2;
		case 4:
			return (phase_fraction(i, j, k) + phase_fraction(i, j, k + 1)) / 2;
		case 5:
			return (phase_fraction(i, j, k) + phase_fraction(i, j, k - 1)) / 2;
	}
}

double effective_viscosity_on_face(int i, int j, int k, int s)
{
	double phase_fraction_on_face_tmp = phase_fraction_on_face(i, j, k, s);
	double shear_rate_on_face_tmp = shear_rate_on_face(i, j, k, s);
	double x = k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= shear_rate_0) {
		x += limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (k_viscosity_snow * pow(shear_rate_on_face_tmp, flow_index - 1) + yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

double pressure_on_face(int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (pressure(i, j, k) + pressure(i + 1, j, k)) / 2;
		case 1:
			return (pressure(i, j, k) + pressure(i - 1, j, k)) / 2;
		case 2:
			return (pressure(i, j, k) + pressure(i, j + 1, k)) / 2;
		case 3:
			return (pressure(i, j, k) + pressure(i, j - 1, k)) / 2;
		case 4:
			return (pressure(i, j, k) + pressure(i, j, k + 1)) / 2;
		case 5:
			return (pressure(i, j, k) + pressure(i, j, k - 1)) / 2;
	}
}

