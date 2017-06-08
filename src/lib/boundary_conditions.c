#include "init_data.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int SET_boundary_CONDITION_velocity_zero_gradient_on_up_and_sides_no_slip_on_low(in I)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * I.stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I.stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	for (p = 0; p < 3; p++) {
		for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
			for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
				for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
					if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
						((!((i >= 0) && (i < I.nx) &&
						   (j >= 0) && (j < I.ny) && 
						   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
						 (k < 0) || (k >= I.nz))) {
							if (k < 0) {
								I.B_prev[B_IND(I, p, i, j, k)] = 0;
							} else if (k >= I.nz) {
								fl_tmp = 0;
								for (l = 0; l < 2 * I.stencil_size + 1; l++) {
									for (m = 0; m < 2 * I.stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < I.nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < I.ny) &&
											(I.ind_cell_multipl[(i + search[l]) * I.ny + j + search[m]] != -1)) {
												I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i + search[l], j + search[m], I.nz - 1)];
												fl_tmp = 1;
												break;
										}
									}
									if (fl_tmp) break;
								}
							} else {
								fl_tmp = 0;
								for (l = 0; l < 2 * I.stencil_size + 1; l++) {
									for (m = 0; m < 2 * I.stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < I.nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < I.ny) &&
											(I.ind_cell_multipl[(i + search[l]) * I.ny + j + search[m]] != -1)) {
												I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i + search[l], j + search[m], k)];
												fl_tmp = 1;
												break;
										}
									}
									if (fl_tmp) break;
								}
							}
					}
				}
			}
		}
	}
	for (p = 0; p < 3; p++) {
		for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
			for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
				for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
					if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
						((!((i >= 0) && (i < I.nx) &&
						   (j >= 0) && (j < I.ny) && 
						   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
						 (k < 0) || (k >= I.nz))) {
							if (k < 0)
								I.B_prev[B_IND(I, p, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_velocity_fixed_value_on_all(in I)
{
	printf("Set the boundary condition for velosity with fixed value on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	for (p = 0; p < 3; p++) {
		for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
			for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
				for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
					if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
						((!((i >= 0) && (i < I.nx) &&
						   (j >= 0) && (j < I.ny) && 
						   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
						 (k < 0) || (k >= I.nz))) {
								I.B_prev[B_IND(I, p, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_phase_fraction_zero_gradient_on_all(in I)
{
	printf("Set the boundary condition for volume snow fraction with zero gradient on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * I.stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I.stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	p = 3;
	for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
		for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
			for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
				if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
					((!((i >= 0) && (i < I.nx) &&
					   (j >= 0) && (j < I.ny) && 
					   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
					 (k < 0) || (k >= I.nz))) {
						if (k < 0) {
							fl_tmp = 0;
							for (l = 0; l < 2 * I.stencil_size + 1; l++) {
								for (m = 0; m < 2 * I.stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < I.nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < I.ny) &&
										(I.ind_cell_multipl[(i + search[l]) * I.ny + j + search[m]] != -1)) {
											I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i + search[l], j + search[m], 0)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						} else if (k >= I.nz) {
							fl_tmp = 0;
							for (l = 0; l < 2 * I.stencil_size + 1; l++) {
								for (m = 0; m < 2 * I.stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < I.nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < I.ny) &&
										(I.ind_cell_multipl[(i + search[l]) * I.ny + j + search[m]] != -1)) {
											I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i + search[l], j + search[m], I.nz - 1)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						} else {
							fl_tmp = 0;
							for (l = 0; l < 2 * I.stencil_size + 1; l++) {
								for (m = 0; m < 2 * I.stencil_size + 1; m++) {
									if ((i + search[l] >= 0) &&
										(i + search[l] < I.nx) &&
										(j + search[m] >= 0) &&
										(j + search[m] < I.ny) &&
										(I.ind_cell_multipl[(i + search[l]) * I.ny + j + search[m]] != -1)) {
											I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i + search[l], j + search[m], k)];
											fl_tmp = 1;
											break;
									}
								}
								if (fl_tmp) break;
							}
						}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_phase_fraction_fixed_value_on_all(in I)
{
	printf("Set the boundary condition for volume snow fraction with fixed values on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	p = 3;
	for (k = -I.stencil_size; k < I.nz + I.stencil_size; k++) {
		for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
			for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
				if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
					((!((i >= 0) && (i < I.nx) &&
					   (j >= 0) && (j < I.ny) && 
					   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
					 (k < 0) || (k >= I.nz))) {
					I.B_prev[B_IND(I, p, i, j, k)] = 0;
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_pressure_fixed_value_on_up_and_sides_zero_gradient_on_low(in I)
{
	printf("Set the boundary condition for pressure with fixed value of atmospheric pressure on upper and sides boundaries, zero gradient condition for lower boundary\n");
	int i, j, k, p;
	p = 4;
	for (k = I.nz + I.stencil_size - 1; k >= - I.stencil_size; k--) {
		for (i = -I.stencil_size; i < I.nx + I.stencil_size; i++) {
			for (j = -I.stencil_size; j < I.ny + I.stencil_size; j++) {
				if ((I.ind_boundary_cells[(i + I.stencil_size) * (I.ny + 2 * I.stencil_size) + j + I.stencil_size] != -1) &&
					((!((i >= 0) && (i < I.nx) &&
					   (j >= 0) && (j < I.ny) && 
					   (I.ind_cell_multipl[i * I.ny + j] != -1))) ||
					 (k < 0) || (k >= I.nz))) {
						if (k >= 0)
							I.B_prev[B_IND(I, p, i, j, k)] = I.pressure_atmosphere;
						else
							I.B_prev[B_IND(I, p, i, j, k)] = I.B_prev[B_IND(I, p, i, j, 0)];
				}
			}
		}
	}
	return 0;
}

