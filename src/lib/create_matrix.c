#include "init_data.h"
#include "utils.h"
#include "x_crank_nikolson_second_combined_VOF.h"
#include "t_second_combined_VOF.h"
#include "create_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DDT(p, i, j, k, object, approximation_order, solution_mode, method) DDT_##object##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define DIV(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) DIV_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define GRAD(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) GRAD_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define VECT(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) VECT_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)

int create_Ab(in I)
{
	printf("Creating matrix A and vector B function\n");
	int i, j, k, p;
	if (I.flag_first_time_step) {
		if ((I.B = (double *) malloc(I.system_dimension * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I.B, 0, I.system_dimension * sizeof(double));
		int num_el_in_row = 100; // it is inaccurate value
		I.non_zero_elem = num_el_in_row * I.system_dimension;
		if ((I.Aelem_csr = (double *) malloc(I.non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I.Aelem_csr, 0, I.non_zero_elem * sizeof(double));
		if ((I.Ajptr_csr = (int *) malloc(I.non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I.Ajptr_csr, -1, I.non_zero_elem * sizeof(int));
		if ((I.Aiptr_csr = (int *) malloc((I.system_dimension + 1) * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I.Aiptr_csr, -1, (I.system_dimension + 1) * sizeof(int));
	} else {
		memset((void *) I.B, 0, I.system_dimension * sizeof(double));
		memset((void *) I.Aelem_csr, 0, I.non_zero_elem * sizeof(double));
	}

/* creating matrix */
	if (I.flag_first_time_step)
		I.A_ind_current = 0;
	for (k = 0; k < I.nz; k++) {
		for (i = 0; i < I.nx; i++) {
			for (j = 0; j < I.ny; j++) {
				if (I.ind_cell_multipl[i * I.ny + j] != -1) {
					/* momentum equation */
					for (p = 0; p < 3; p++) {
						if (I.flag_first_time_step) {
							if (I.Ajptr_csr[I.A_ind_current] != -1)
								I.A_ind_current++;
							I.Aiptr_csr[A_IND(I, p, i, j, k)] = I.A_ind_current;
						}
						if (DDT(p, i, j, k, density_velocity, second, combined, VOF)) return 1;
						//if (DDT(p, i, j, k, density_velocity)) return 1;
						if (DIV(p, i, j, k, density_velocity_velocity, crank_nikolson, second, combined, VOF)) return 1;
						//if (VECT(p, i, j, k, gravity_force, crank_nikolson, second, combined, VOF)) return 1;
						if (GRAD(p, i, j, k, pressure, crank_nikolson, second, combined, VOF)) return 1;
						if (DIV(p, i, j, k, shear_stress, crank_nikolson, second, combined, VOF)) return 1;
					}
					/* transport equation for snow volume fraction */
					p = 3;
					if (I.flag_first_time_step) {
						if (I.Ajptr_csr[I.A_ind_current] != -1)
							I.A_ind_current++;
						I.Aiptr_csr[A_IND(I, p, i, j, k)] = I.A_ind_current;
					}
					if (DDT(p, i, j, k, snow_volume_fraction, second, combined, VOF)) return 1;
					if (DIV(p, i, j, k, snow_volume_fraction_velocity, crank_nikolson, second, combined, VOF)) return 1;
					//if (DIV(p, i, j, k, grad_snow_volume_fraction, crank_nikolson, second, combined, FDM)) return 1;
					/* continuity equation */
					p = 4;
					if (I.flag_first_time_step) {
						if (I.Ajptr_csr[I.A_ind_current] != -1)
							I.A_ind_current++;
						I.Aiptr_csr[A_IND(I, p, i, j, k)] = I.A_ind_current;
					}
					if (DIV(p, i, j, k, velocity_cont, crank_nikolson, second, combined, VOF)) return 1;
					/* poisson equation for pressure */
					//p = 4;
					//if (flag_first_time_step) {
					//	if (Ajptr_csr[A_ind_current] != -1)
					//		A_ind_current++;
					//	Aiptr_csr[A_IND(p, i, j, k)] = A_ind_current;
					//}
					//if (DIV(p, i, j, k, grad_pressure, crank_nikolson, second, combined, VOF)) return 1;
					//if (DIV(p, i, j, k, div_density_velocity_velocity, crank_nikolson, second, combined, VOF)) return 1;
				}
			}
		}
	}
	if (I.flag_first_time_step) {
		I.non_zero_elem = I.A_ind_current + 1;
		I.Aiptr_csr[I.system_dimension] = I.non_zero_elem;
		if ((I.Aelem_csr = (double *) realloc(I.Aelem_csr, I.non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		if ((I.Ajptr_csr = (int *) realloc(I.Ajptr_csr, I.non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
	}
	//if (flag_first_time_step)
	//	print_A_csr();
	return 0;
}

