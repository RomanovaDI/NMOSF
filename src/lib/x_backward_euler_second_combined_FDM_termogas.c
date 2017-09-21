#include "init_data.h"
#include "utils.h"
#include "x_backward_euler_second_combined_FDM_termogas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}
#define DDT(p, i, j, k, object, approximation_order, solution_mode, method, task) DDT_##object##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define DIV(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) DIV_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define GRAD(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) GRAD_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define VECT(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) VECT_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define SCAL(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) SCAL_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)
#define LAPL(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method, task) LAPL_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method##_##task(I, p, i, j, k)

int DIV_concentration_density_average_velocity_backward_euler_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3];
	double A_value;
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		A_value = density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (2 * I->dx[pr]);
		WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		A_value = - density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
		WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
	}
	return 0;
}

int SCAL_mass_inflow_rate_backward_euler_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k);
	return 0;
}

int LAPL_coef_pressure_backward_euler_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		A_value = - Darsi_M_coef / (I->dx[pr] * I->dx[pr]);
		WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		A_value *= -2;
		WRITE_TO_A(p, i, j, k, -1);
	}
	return 0;
}

int DIV_density_average_velocity_backward_euler_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		I->B[A_IND(I, p, i, j, k)] -= (density_t(I, p - 4, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, p - 4, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
			density_t(I, p - 4, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, p - 4, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
	}
	return 0;
}

int DIV_density_saturation_internal_energy_avarage_velocity_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pp, pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		A_value = 0;
		for (pp = 0; pp < 2; pp++)
			A_value += density_t(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				saturation(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				I->specific_heat[pp] *
				avarage_velocity(I, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
		for (pp = 0; pp < 4; pp++)
			A_value += density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				saturation(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				I->specific_heat[pp + 2] *
				concentration(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
		A_value /= 2 * I->dx[pr];
		WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		A_value = 0;
		for (pp = 0; pp < 2; pp++)
			A_value += density_t(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				saturation(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				I->specific_heat[pp] *
				avarage_velocity(I, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
		for (pp = 0; pp < 4; pp++)
			A_value += density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				saturation(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				I->specific_heat[pp + 2] *
				concentration(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
		A_value /= 2 * I->dx[pr];
		WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
	}
}

int LAPL_heat_influx_vector_backward_euler_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = I->porousness * 
	return 0;
}
