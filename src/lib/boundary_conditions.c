#include "init_data.h"
#include "utils.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include "boundary_conditions.h"

#if AVALANCHE
int SET_boundary_CONDITION_velocity_zero_gradient_on_up_and_sides_no_slip_on_low(in *I)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * I->stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	for (p = 0; p < 3; p++) {
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
						((!((i >= 0) && (i < I->nx) &&
						   (j >= 0) && (j < I->ny) && 
						   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
						 (k < 0) || (k >= I->nz))) {
							if (k < 0) {
								I->B_prev[B_IND(I, p, i, j, k)] = 0;
							} else if (k >= I->nz) {
								fl_tmp = 0;
								for (l = 0; l < 2 * I->stencil_size + 1; l++) {
									for (m = 0; m < 2 * I->stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < I->nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < I->ny) &&
											(I->ind_cell_multipl[(i + search[l]) * I->ny + j + search[m]] != -1)) {
												I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i + search[l], j + search[m], I->nz - 1)];
												fl_tmp = 1;
												break;
										}
									}
									if (fl_tmp) break;
								}
							} else {
								fl_tmp = 0;
								for (l = 0; l < 2 * I->stencil_size + 1; l++) {
									for (m = 0; m < 2 * I->stencil_size + 1; m++) {
										if ((i + search[l] >= 0) &&
											(i + search[l] < I->nx) &&
											(j + search[m] >= 0) &&
											(j + search[m] < I->ny) &&
											(I->ind_cell_multipl[(i + search[l]) * I->ny + j + search[m]] != -1)) {
												I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i + search[l], j + search[m], k)];
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
	//for (p = 0; p < 3; p++) {
	//	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
	//		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
	//			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
	//				if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
	//					((!((i >= 0) && (i < I->nx) &&
	//					   (j >= 0) && (j < I->ny) && 
	//					   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
	//					 (k < 0) || (k >= I->nz))) {
	//						if (k < 0)
	//							I->B_prev[B_IND(I, p, i, j, k)] = 0;
	//				}
	//			}
	//		}
	//	}
	//}
	return 0;
}

int SET_boundary_CONDITION_velocity_zero_gradient_on_y_sides_no_slip_on_other_upper_wall_is_mooving(in *I)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, m;
	int search[2 * I->stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	for (p = 0; p < 3; p++) {
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if (boundary_cell(I, i, j, k)) {
						if ((j < 0) || (j >= I->ny)) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((j + search[m] >= 0) &&
									(j + search[m] < I->ny) &&
									(I->ind_cell_multipl[i * I->ny + j + search[m]] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j + search[m], k)];
										break;
								}
							}
						} else if ((k >= I->nz) && (p == 0)) {
							I->B_prev[B_IND(I, p, i, j, k)] = 0.1;
						} else {
							I->B_prev[B_IND(I, p, i, j, k)] = 0;
						}
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_velocity_zero_gradient_on_y_and_x_sides_no_slip_on_other_upper_wall_is_mooving(in *I)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, m;
	int search[2 * I->stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	for (p = 0; p < 3; p++) {
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if (boundary_cell(I, i, j, k)) {
						//if (((j < 0) || (j >= I->ny)) && (i >= 0) && (i < I->nx)) {
						if ((j < 0) || (j >= I->ny)) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((j + search[m] >= 0) &&
									(j + search[m] < I->ny) &&
									(I->ind_cell_multipl[i * I->ny + j + search[m]] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j + search[m], k)];
										break;
								}
							}
						//} else if (((i < 0) || (i >= I->nx)) && (j >= 0) && (j < I->ny)) {
						} else if ((i < 0) || (i >= I->nx)) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((i + search[m] >= 0) &&
									(i + search[m] < I->nx) &&
									(I->ind_cell_multipl[(i + search[m]) * I->ny + j] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i + search[m], j, k)];
										break;
								}
							}
						} else if ((k >= I->nz) && (p == 0)) {
							I->B_prev[B_IND(I, p, i, j, k)] = 1;
						} else {
							I->B_prev[B_IND(I, p, i, j, k)] = 0;
						}
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_velocity_zero_gradient_on_y_and_x_and_upper_sides_no_slip_on_low(in *I)
{
	printf("Set the boundary condition for velosity with zero gradient on upper and sides boundaries and no slip condition on the lower boundary\n");
	int i, j, k, p, m;
	int search[2 * I->stencil_size + 1];
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	for (p = 0; p < 3; p++) {
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if (boundary_cell(I, i, j, k)) {
						if ((j < 0) || (j >= I->ny)) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((j + search[m] >= 0) &&
									(j + search[m] < I->ny) &&
									(I->ind_cell_multipl[i * I->ny + j + search[m]] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j + search[m], k)];
										break;
								}
							}
						} else if ((i < 0) || (i >= I->nx)) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((i + search[m] >= 0) &&
									(i + search[m] < I->nx) &&
									(I->ind_cell_multipl[(i + search[m]) * I->ny + j] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i + search[m], j, k)];
										break;
								}
							}
						} else if (k >= I->nz) {
							for (m = 0; m < 2 * I->stencil_size + 1; m++) {
								if ((k + search[m] >= 0) &&
									(k + search[m] < I->nz) &&
									(I->ind_cell_multipl[i * I->ny + j] != -1)) {
										I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, k + search[m])];
										break;
								}
							}
						} else {
							I->B_prev[B_IND(I, p, i, j, k)] = 0;
						}
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_velocity_fixed_value_on_all(in *I)
{
	printf("Set the boundary condition for velosity with fixed value on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	for (p = 0; p < 3; p++) {
		for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
						((!((i >= 0) && (i < I->nx) &&
						   (j >= 0) && (j < I->ny) && 
						   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
						 (k < 0) || (k >= I->nz))) {
								I->B_prev[B_IND(I, p, i, j, k)] = 0;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_phase_fraction_zero_gradient_on_all(in *I)
{
	printf("Set the boundary condition for volume snow fraction with zero gradient on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * I->stencil_size + 1];
	double x;
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	p = 3;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (boundary_cell(I, i, j, k)) {
					if (k < 0)
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, 0)];
					else if (k >= I->nz)
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, I->nz - 1)];
					else if (count_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 1, j, k)) && (! boundary_cell(I, i - 1, j, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j, k)) && (! boundary_cell(I, i + 1, j, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j - 1, k)) && (! boundary_cell(I, i, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j + 1, k)) && (! boundary_cell(I, i, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i, j + 1, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					} else if (count_second_order_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 2, j, k)) && (! boundary_cell(I, i - 2, j, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j, k)) && (! boundary_cell(I, i + 2, j, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j - 2, k)) && (! boundary_cell(I, i, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j + 2, k)) && (! boundary_cell(I, i, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i, j + 2, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					} else if (count_other_corner_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 1, j - 1, k)) && (! boundary_cell(I, i - 1, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j + 1, k)) && (! boundary_cell(I, i - 1, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j - 1, k)) && (! boundary_cell(I, i + 1, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j + 1, k)) && (! boundary_cell(I, i + 1, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j - 2, k)) && (! boundary_cell(I, i - 1, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j + 2, k)) && (! boundary_cell(I, i - 1, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j - 2, k)) && (! boundary_cell(I, i + 1, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j + 2, k)) && (! boundary_cell(I, i + 1, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j - 1, k)) && (! boundary_cell(I, i - 2, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j + 1, k)) && (! boundary_cell(I, i - 2, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j - 1, k)) && (! boundary_cell(I, i + 2, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j + 1, k)) && (! boundary_cell(I, i + 2, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j - 2, k)) && (! boundary_cell(I, i - 2, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j + 2, k)) && (! boundary_cell(I, i - 2, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j - 2, k)) && (! boundary_cell(I, i + 2, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j + 2, k)) && (! boundary_cell(I, i + 2, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j + 2, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					}
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_phase_fraction_fixed_value_on_all(in *I)
{
	printf("Set the boundary condition for volume snow fraction with fixed values on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	p = 3;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
					((!((i >= 0) && (i < I->nx) &&
					   (j >= 0) && (j < I->ny) && 
					   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
					 (k < 0) || (k >= I->nz))) {
					I->B_prev[B_IND(I, p, i, j, k)] = 0;
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_pressure_fixed_value_on_up_and_sides_zero_gradient_on_low(in *I)
{
	printf("Set the boundary condition for pressure with fixed value of atmospheric pressure on upper and sides boundaries, zero gradient condition for lower boundary\n");
	int i, j, k, p;
	p = 4;
	for (k = I->nz + I->stencil_size - 1; k >= - I->stencil_size; k--) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
					((!((i >= 0) && (i < I->nx) &&
					   (j >= 0) && (j < I->ny) && 
					   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
					 (k < 0) || (k >= I->nz))) {
						if (k >= 0)
							I->B_prev[B_IND(I, p, i, j, k)] = I->pressure_atmosphere;
						else
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, 0)];
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_pressure_zero_gradient_on_all(in *I)
{
	printf("Set the boundary condition for volume snow fraction with zero gradient on all boundaries\n");
	int i, j, k, p, l, m, fl_tmp;
	int search[2 * I->stencil_size + 1];
	double x;
	search[0] = 0;
	for (i = 1; i <= I->stencil_size; i++) {
		search[i * 2 - 1] = i;
		search[i * 2] = -i;
	}
	p = 4;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (boundary_cell(I, i, j, k)) {
					if (k < 0)
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, 0)];
					else if (k >= I->nz)
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, I->nz - 1)];
					else if (count_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 1, j, k)) && (! boundary_cell(I, i - 1, j, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j, k)) && (! boundary_cell(I, i + 1, j, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j - 1, k)) && (! boundary_cell(I, i, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j + 1, k)) && (! boundary_cell(I, i, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i, j + 1, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					} else if (count_second_order_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 2, j, k)) && (! boundary_cell(I, i - 2, j, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j, k)) && (! boundary_cell(I, i + 2, j, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j - 2, k)) && (! boundary_cell(I, i, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i, j + 2, k)) && (! boundary_cell(I, i, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i, j + 2, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					} else if (count_other_corner_neighbor_internal_cells(I, i, j, k)) {
						x = l = 0;
						if ((cell_of_computation_domain(I, i - 1, j - 1, k)) && (! boundary_cell(I, i - 1, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j + 1, k)) && (! boundary_cell(I, i - 1, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j - 1, k)) && (! boundary_cell(I, i + 1, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j + 1, k)) && (! boundary_cell(I, i + 1, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j - 2, k)) && (! boundary_cell(I, i - 1, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 1, j + 2, k)) && (! boundary_cell(I, i - 1, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 1, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j - 2, k)) && (! boundary_cell(I, i + 1, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 1, j + 2, k)) && (! boundary_cell(I, i + 1, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 1, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j - 1, k)) && (! boundary_cell(I, i - 2, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j + 1, k)) && (! boundary_cell(I, i - 2, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j - 1, k)) && (! boundary_cell(I, i + 2, j - 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j - 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j + 1, k)) && (! boundary_cell(I, i + 2, j + 1, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j + 1, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j - 2, k)) && (! boundary_cell(I, i - 2, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i - 2, j + 2, k)) && (! boundary_cell(I, i - 2, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i - 2, j + 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j - 2, k)) && (! boundary_cell(I, i + 2, j - 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j - 2, k)];
							l++;
						}
						if ((cell_of_computation_domain(I, i + 2, j + 2, k)) && (! boundary_cell(I, i + 2, j + 2, k))) {
							x += I->B_prev[B_IND(I, p, i + 2, j + 2, k)];
							l++;
						}
						x /= (double) l;
						I->B_prev[B_IND(I, p, i, j, k)] = x;
					}
				}
			}
		}
	}
	return 0;
}
#endif

#if TERMOGAS
int SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out_consistent(in *I)
{
#if DEBUG
	printf("Set the boundary condition in termogas case for all values with no boundaries condition, 4 injection well and 1 production well.\n");
#endif
	int i, j, k, p, l, m, fl_tmp;
	for (i = 0; i < I->nx; i++) {
		for (j = 0; j < I->ny; j++) {
			if (injection_well(I, i, j, 0))
				set_injection_well(I, i, j, 0);
			if (production_well(I, i, j, 0))
				set_production_well(I, i, j, 0);
		}
	}
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (boundary_cell(I, i, j, k)) {
					for (p = 0; p < I->num_parameters; p++) {
						//printf("B_IND = %d, B_prev = %lf; INT_NEI_IND = %d, B_prev = %lf\n", B_IND(I, p, i, j, k), I->B_prev[B_IND(I, p, i, j, k)], INT_NEI_IND(I, p, i, j, k), I->B_prev[INT_NEI_IND(I, p, i, j, k)]);
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[INT_NEI_IND(I, p, i, j, k)];
					}
				}
			}
		}
	}
	/*
	for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
		for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
			if (boundary_cell(I, i, j, 0)) {
				for (p = 0; p < I->num_parameters; p++) {
					if (i < 0) {
						if (j < 0)
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, 0, 0, 0)];
						else if (j >= I->ny)
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, 0, I->ny - 1, 0)];
						else
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, 0, j, 0)];
					} else if (i >= I->nx) {
						if (j < 0)
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, I->nx - 1, 0, 0)];
						else if (j >= I->ny)
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, I->nx - 1, I->ny - 1, 0)];
						else
							I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, I->nx - 1, j, 0)];
					} else if (j < 0)
						I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, i, 0, 0)];
					else if (j >= I->ny)
						I->B_prev[B_IND(I, p, i, j, 0)] = I->B_prev[B_IND(I, p, i, I->ny - 1, 0)];
					I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[INT_NEI_IND(I, p, i, j, k)];
				}
			}
		}
	}
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		if (k != 0) {
			for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
				for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
					for (p = 0; p < I->num_parameters; p++) {
						I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, 0)];
					}
				}
			}
		}
	}
	*/
	return 0;
}

int set_injection_well(in *I, int i, int j, int k)
{
	/*
	for (int ii = -1; ii < 2; ii++) {
		for (int jj = -1; jj < 2; jj++) {
			for (int kk = -1; kk < 2; kk++) {
				I->B_prev[B_IND(I, 0, i + ii, j + jj, k + kk)] = 3.55 * (1 - 2 * I->epsilon) / 4.55;
				I->B_prev[B_IND(I, 1, i + ii, j + jj, k + kk)] = (1 - 2 * I->epsilon) / 4.55;
				I->B_prev[B_IND(I, 2, i + ii, j + jj, k + kk)] = I->epsilon;
				I->B_prev[B_IND(I, 3, i + ii, j + jj, k + kk)] = I->epsilon;
				I->B_prev[B_IND(I, 5, i + ii, j + jj, k + kk)] = 0.5 - (I->residual_saturation[1] + I->epsilon) / 2;
				I->B_prev[B_IND(I, 6, i + ii, j + jj, k + kk)] = I->residual_saturation[1] + I->epsilon;
				I->B_prev[B_IND(I, 7, i + ii, j + jj, k + kk)] = 0.5 - (I->residual_saturation[1] + I->epsilon) / 2;
				I->B_prev[B_IND(I, 8, i + ii, j + jj, k + kk)] = I->injection_well_temperature;
			}
		}
	}
	*/
	I->B_prev[B_IND(I, 0, i, j, k)] = 3.55 / 4.55;
	//I->B_prev[B_IND(I, 0, i, j, k)] = 3.55 * (1 - 2 * I->epsilon) / 4.55;
	I->B_prev[B_IND(I, 1, i, j, k)] = 1 / 4.55;
	//I->B_prev[B_IND(I, 1, i, j, k)] = (1 - 2 * I->epsilon) / 4.55;
	I->B_prev[B_IND(I, 2, i, j, k)] = 0;
	I->B_prev[B_IND(I, 3, i, j, k)] = 0;
	I->B_prev[B_IND(I, 4, i, j, k)] = I->injection_well_pressure;
	I->B_prev[B_IND(I, 5, i, j, k)] = 0.7 - I->residual_saturation[1] / 2 - I->residual_saturation[0] - 0.05;
	I->B_prev[B_IND(I, 6, i, j, k)] = 0.05;
	I->B_prev[B_IND(I, 7, i, j, k)] = 0.3 - I->residual_saturation[1] / 2 - I->residual_saturation[2];
	//I->B_prev[B_IND(I, 7, i, j, k)] = I->residual_saturation[2] + I->epsilon;//0.3 - (I->residual_saturation[1] + I->epsilon) / 2;
	I->B_prev[B_IND(I, 8, i, j, k)] = I->injection_well_temperature - I->initial_temperature;
	return 0;
}

int set_production_well(in *I, int i, int j, int k)
{
	I->B_prev[B_IND(I, 4, i, j, k)] = I->production_well_pressure;
	//I->B_prev[B_IND(I, 8, i, j, k)] = I->initial_temperature;
	//I->B_prev[B_IND(I, 4, i+1, j, k)] = I->production_well_pressure;
	//I->B_prev[B_IND(I, 4, i-1, j, k)] = I->production_well_pressure;
	return 0;
}

int SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out_parallel(in *I)
{
#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Set the boundary condition in termogas case for all values with no boundaries condition, 4 injection well and 1 production well in process %d.\n", I->my_rank);
#endif
	//reconstruct_src(I);
	MPI_Bcast(I->gl_B, I->gl_n_cells_multipl * I->gl_nz * I->num_parameters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int i, j, k, p, l, m, fl_tmp;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if ((internal_cell(I, i, j, k)) && (I->ind_proc[(i + I->ind_start_region_proc[0]) * I->gl_ny + j + I->ind_start_region_proc[1]] != I->my_rank)) {
					for (p = 0; p < I->num_parameters; p++) {
						I->B_prev[B_IND(I, p, i, j, k)] = I->gl_B[GL_A_IND(I, p, i + I->ind_start_region_proc[0], j + I->ind_start_region_proc[1], k)];
					}
				} else if (boundary_cell(I, i, j, k)) {
					for (p = 0; p < I->num_parameters; p++) {
						if (k < 0)
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, 0)];
						else if (k >= I->nz)
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, I->nz - 1)];
						else if ((i < 0) && (I->my_rank < I->y_regions))
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, 0, j, k)];
						else if ((i < 0) && (I->my_rank >= I->y_regions))
							I->B_prev[B_IND(I, p, i, j, k)] = I->gl_B[GL_A_IND(I, p, i + I->ind_start_region_proc[0], j + I->ind_start_region_proc[1], k)];
						else if ((i >= I->nx) && (I->my_rank >= I->y_regions * (I->x_regions - 1)))
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, I->nx - 1, j, k)];
						else if ((i >= I->nx) && (I->my_rank < I->y_regions * (I->x_regions - 1)))
							I->B_prev[B_IND(I, p, i, j, k)] = I->gl_B[GL_A_IND(I, p, i + I->ind_start_region_proc[0], j + I->ind_start_region_proc[1], k)];
						else if ((j < 0) && (I->my_rank % I->y_regions == 0))
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, 0, k)];
						else if ((j < 0) && (I->my_rank % I->y_regions != 0))
							I->B_prev[B_IND(I, p, i, j, k)] = I->gl_B[GL_A_IND(I, p, i + I->ind_start_region_proc[0], j + I->ind_start_region_proc[1], k)];
						else if ((j >= I->ny) && (I->my_rank % I->y_regions == I->y_regions - 1))
							I->B_prev[B_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, I->ny - 1, k)];
						else if ((j >= I->ny) && (I->my_rank % I->y_regions != I->y_regions - 1))
							I->B_prev[B_IND(I, p, i, j, k)] = I->gl_B[GL_A_IND(I, p, i + I->ind_start_region_proc[0], j + I->ind_start_region_proc[1], k)];
					}
				}
			}
		}
	}
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (injection_well(I, i, j, k))
					set_injection_well(I, i, j, k);
				if (production_well(I, i, j, k)) {
					set_production_well(I, i, j, k);
				}
			}
		}
	}
	return 0;
}

int SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out(in *I)
{
	if (I->nproc > 1) {
		if (SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out_parallel(I)) return 1;
	} else {
		if (SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out_consistent(I)) return 1;
	}
	return 0;
}
#endif
