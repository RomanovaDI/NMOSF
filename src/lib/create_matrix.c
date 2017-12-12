#include "init_data.h"
#include "utils.h"
#include "matrix_functions.h"
#include "mpi.h"
#if AVALANCHE
#include "x_crank_nikolson_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_FDM_avalanche.h"
#include "x_backward_euler_second_combined_VOF_avalanche.h"
#include "x_backward_euler_second_ultra_combined_VOF_avalanche.h"
#include "t_second_combined_VOF_avalanche.h"
#include "t_second_ultra_combined_VOF_avalanche.h"
#endif
#if TERMOGAS
#include "t_second_combined_FDM_termogas.h"
#include "x_backward_euler_second_combined_FDM_termogas.h"
#include "t_second_separated_FDM_termogas.h"
#include "x_backward_euler_second_separated_FDM_termogas.h"
#endif
#include "t_test.h"
#include "create_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#define DDT(p, i, j, k, object, approximation_order, solution_mode, method, task) DDT_##object##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define DIV(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) DIV_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define GRAD(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) GRAD_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define VECT(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) VECT_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define SCAL(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) SCAL_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define LAPL(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) LAPL_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)

int create_Ab_avalanche(in *I)
{
#if AVALANCHE
	printf("Creating matrix A and vector B function\n");
	int i, j, k, p;
	if (I->flag_first_time_step) {
		if ((I->B = (double *) malloc(I->system_dimension * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->B, 0, I->system_dimension * sizeof(double));
		int num_el_in_row = 100; // it is inaccurate value
		I->non_zero_elem = num_el_in_row * I->system_dimension;
		if ((I->Aelem_csr = (double *) malloc(I->non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Aelem_csr, 0, I->non_zero_elem * sizeof(double));
		if ((I->Ajptr_csr = (int *) malloc(I->non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Ajptr_csr, -1, I->non_zero_elem * sizeof(int));
		if ((I->Aiptr_csr = (int *) malloc((I->system_dimension + 1) * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Aiptr_csr, -1, (I->system_dimension + 1) * sizeof(int));
	} else {
		memset((void *) I->B, 0, I->system_dimension * sizeof(double));
		memset((void *) I->Aelem_csr, 0, I->non_zero_elem * sizeof(double));
	}

/* creating matrix */
	if (I->flag_first_time_step)
		I->A_ind_current = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					/* momentum equation */
					for (p = 0; p < 3; p++) {
						if (I->flag_first_time_step) {
							if (I->Ajptr_csr[I->A_ind_current] != -1)
								I->A_ind_current++;
							I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
						}
						if (DDT(p, i, j, k, density_velocity, second, combined, VOF, avalanche)) return 1;
						if (DIV(p, i, j, k, density_velocity_velocity, crank_nikolson, second, combined, VOF, avalanche)) return 1;
						//if (VECT(p, i, j, k, gravity_force, backward_euler, second, combined, VOF, avalanche)) return 1;
						//if (GRAD(p, i, j, k, pressure, crank_nikolson, second, combined, VOF, avalanche)) return 1;
						if (GRAD(p, i, j, k, pressure, crank_nikolson, second, combined, VOF, avalanche)) return 1;
						if (DIV(p, i, j, k, shear_stress_linear, crank_nikolson, second, combined, VOF, avalanche)) return 1;
						//if (DIV(p, i, j, k, shear_stress, crank_nikolson, second, combined, VOF, avalanche)) return 1;
					}
					/* transport equation for snow volume fraction */
					p = 3;
					if (I->flag_first_time_step) {
						if (I->Ajptr_csr[I->A_ind_current] != -1)
							I->A_ind_current++;
						I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					}
					//if (VECT(p, i, j, k, barotropy_density, forward_euler, second, combined, VOF, avalanche)) return 1;
					if (DDT(p, i, j, k, density, second, combined, VOF, avalanche)) return 1;
					if (DIV(p, i, j, k, density_velocity, crank_nikolson, second, combined, VOF, avalanche)) return 1;
					//if (DIV(p, i, j, k, grad_snow_volume_fraction, crank_nikolson, second, combined, FDM, avalanche)) return 1;
					/* continuity equation */
					//p = 4;
					//if (I->flag_first_time_step) {
					//	if (I->Ajptr_csr[I->A_ind_current] != -1)
					//		I->A_ind_current++;
					//	I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					//}
					//if (DIV(p, i, j, k, velocity_cont, crank_nikolson, second, combined, VOF, avalanche)) return 1;
					/* poisson equation for pressure */
					p = 4;
					if (I->flag_first_time_step) {
						if (I->Ajptr_csr[I->A_ind_current] != -1)
							I->A_ind_current++;
						I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					}
					////if (DDT(p, i, j, k, pressure_cont, second, combined, VOF, avalanche)) return 1;
					////if (DIV(p, i, j, k, snow_volume_fraction_velocity, forward_euler, second, combined, VOF, avalanche)) return 1;
					if (VECT(p, i, j, k, barotropy_pressure, crank_nikolson, second, combined, VOF, avalanche)) return 1;
					//if (DIV(p, i, j, k, grad_pressure, backward_euler, second, combined, FDM, avalanche)) return 1;
					//if (DIV(p, i, j, k, div_density_velocity_velocity, forward_euler, second, combined, FDM, avalanche)) return 1;
				}
			}
		}
	}
	if (I->flag_first_time_step) {
		I->non_zero_elem = I->A_ind_current + 1;
		I->Aiptr_csr[I->system_dimension] = I->non_zero_elem;
		if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
	}
	//if (flag_first_time_step)
	//	print_A_csr();
#endif
	return 0;
}

int create_Ab_termogas(in *I)
{
#if TERMOGAS
#if DEBUG
	printf("Creating matrix A and vector B function in process %d\n", I->my_rank);
#endif
	int i, j, k, p;
	if (I->flag_first_time_step) {
		if ((I->B = (double *) malloc(I->system_dimension * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->B, 0, I->system_dimension * sizeof(double));
		int num_el_in_row = (2 * I->stencil_size * 3 + 1) * I->num_parameters;
		I->non_zero_elem = num_el_in_row * I->system_dimension;
		if ((I->Aelem_csr = (double *) malloc(I->non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Aelem_csr, 0, I->non_zero_elem * sizeof(double));
		if ((I->Ajptr_csr = (int *) malloc(I->non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Ajptr_csr, -1, I->non_zero_elem * sizeof(int));
		if ((I->Aiptr_csr = (int *) malloc((I->system_dimension + 1) * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		memset((void *) I->Aiptr_csr, -1, (I->system_dimension + 1) * sizeof(int));
	} else {
		memset((void *) I->B, 0, I->system_dimension * sizeof(double));
		memset((void *) I->Aelem_csr, 0, I->non_zero_elem * sizeof(double));
	}

/* creating matrix */
	if (I->flag_first_time_step)
		I->A_ind_current = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					/* concentration equation */
					for (p = 0; p < 4; p++) {
						if (I->flag_first_time_step) {
							if (I->Ajptr_csr[I->A_ind_current] != -1)
								I->A_ind_current++;
							I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
						}
						//if (I->equation_num == 0) {
						if (1) {
							if (production_well(I, i, j, k)) {
								if (DDT(p, i, j, k, arithmetic_mean_of_neighboring_cells, second, separated, FDM, termogas)) return 0;
							} else if (injection_well(I, i, j, k)) {
								if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
							} else {
								if (DDT(p, i, j, k, concentration_density_saturation_porousness, second, separated, FDM, termogas)) return 1;
								if (DIV(p, i, j, k, concentration_density_average_velocity, backward_euler, second, separated, FDM, termogas)) return 1;
								if (SCAL(p, i, j, k, mass_inflow_rate, backward_euler, second, separated, FDM, termogas)) return 1;
							}
						} else {
							if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
						}
					}
					/* pressure equation */
					p = 4;
					if (I->flag_first_time_step) {
						if (I->Ajptr_csr[I->A_ind_current] != -1)
							I->A_ind_current++;
						I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					}
					//if ((I->equation_num == 4) && (!(well(I, i, j, k)))) {
					if ((1) && (!(well(I, i, j, k)))) {
						if (DDT(p, i, j, k, coef_pressure, second, separated, FDM, termogas)) return 1;
						if (LAPL(p, i, j, k, coef_pressure, backward_euler, second, separated, FDM, termogas)) return 1;
					} else {
						if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
					}
					/* saturation equation */
					for (p = 5; p < 8; p++) {
						if (I->flag_first_time_step) {
							if (I->Ajptr_csr[I->A_ind_current] != -1)
								I->A_ind_current++;
							I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
						}
						//if (I->equation_num == 1) {
						if (1) {
							if (production_well(I, i, j, k)) {
								if (DDT(p, i, j, k, arithmetic_mean_of_neighboring_cells, second, separated, FDM, termogas)) return 0;
								//if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
							} else if (injection_well(I, i, j, k)) {
								if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
							} else {
								if (DDT(p, i, j, k, density_saturation_porousness, second, separated, FDM, termogas)) return 1;
								if (DIV(p, i, j, k, density_average_velocity, backward_euler, second, separated, FDM, termogas)) return 1;
								if (SCAL(p, i, j, k, mass_inflow_rate, backward_euler, second, separated, FDM, termogas)) return 1;
							}
						} else {
							if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
						}
					}
					/* temperature flow equation */
					p = 8;
					if (I->flag_first_time_step) {
						if (I->Ajptr_csr[I->A_ind_current] != -1)
							I->A_ind_current++;
						I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					}
					//if (I->equation_num == 2) {
					if (1) {
						if (production_well(I, i, j, k)) {
							if (DDT(p, i, j, k, arithmetic_mean_of_neighboring_cells, second, separated, FDM, termogas)) return 0;
						} else if (injection_well(I, i, j, k)) {
							if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
						} else {
							if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
							if (DDT(p, i, j, k, porousness_density_energy_flow, second, separated, FDM, termogas)) return 1;
							if (DIV(p, i, j, k, density_saturation_internal_energy_avarage_velocity, backward_euler, second, separated, FDM, termogas)) return 1;
							if (DIV(p, i, j, k, heat_influx_vector_flow, backward_euler, second, separated, FDM, termogas)) return 1;
							if (SCAL(p, i, j, k, heat_flow, backward_euler, second, separated, FDM, termogas)) return 1;
							if (SCAL(p, i, j, k, chemical_reaction_heat_flow, backward_euler, second, separated, FDM, termogas)) return 1;
						}
					} else {
						if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
					}
					/* temperature environment equation */
					p = 9;
					if (I->flag_first_time_step) {
						if (I->Ajptr_csr[I->A_ind_current] != -1)
							I->A_ind_current++;
						I->Aiptr_csr[A_IND(I, p, i, j, k)] = I->A_ind_current;
					}
					//if (I->equation_num == 3) {
					if (1) {
						if (production_well(I, i, j, k)) {
							if (DDT(p, i, j, k, arithmetic_mean_of_neighboring_cells, second, separated, FDM, termogas)) return 0;
						} else if (injection_well(I, i, j, k)) {
							if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
						} else {
							if (DDT(p, i, j, k, porousness_density_energy_environment, second, separated, FDM, termogas)) return 1;
							if (DIV(p, i, j, k, heat_influx_vector_environment, backward_euler, second, separated, FDM, termogas)) return 1;
							if (SCAL(p, i, j, k, minus_heat_flow, backward_euler, second, separated, FDM, termogas)) return 1;
						}
					} else {
						if (DDT(p, i, j, k, identity, second, separated, FDM, termogas)) return 1;
					}
				}
			}
		}
	}
	if (I->flag_first_time_step) {
		I->non_zero_elem = I->A_ind_current + 1;
		I->Aiptr_csr[I->system_dimension] = I->non_zero_elem;
		if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
		if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
			printf("Memory error\n");
			return 1;
		}
	}
	// (I->equation_num == 1)
	//if (I->flag_first_time_step)
	//	print_A_csr(I);
#endif
	return 0;
}
