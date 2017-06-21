#include "init_data.h"
#include "utils.h"
#include "x_forward_euler_second_combined_FDM.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}
#define DDT(p, i, j, k, object, approximation_order, solution_mode, method) DDT_##object##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define DIV(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) DIV_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define GRAD(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) GRAD_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)
#define VECT(p, i, j, k, object, numerical_scheme, approximation_order, solution_mode, method) VECT_##object##_##numerical_scheme##_##approximation_order##_##solution_mode##_##method(I, p, i, j, k)

int GRAD_pressure_forward_euler_second_combined_FDM(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	I->B[A_IND(I, p, i, j, k)] -= (pressure(I, i + ind_p[0], j + ind_p[1], k + ind_p[2]) - pressure(I, i - ind_p[0], j - ind_p[1], k - ind_p[2])) / (2 * I->dx[p]);
	return 0;
}

int DIV_density_velocity_velocity_forward_euler_second_combined_FDM(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[p] = 1;
		I->B[A_IND(I, p, i, j, k)] -= density(I, i, j, k) * velocity(I, pr, i, j, k) *
			(velocity(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - velocity(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
	}
	return 0;
}

int DIV_snow_volume_fraction_velocity_forward_euler_second_combined_FDM(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[p] = 1;
		I->B[A_IND(I, p, i, j, k)] -= (phase_fraction(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * velocity(I, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
				phase_fraction(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * velocity(I, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
	}
	return 0;
}

int DIV_shear_stress_linear_forward_euler_second_combined_FDM(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3], ind_p[3];
	//for (pr = 0; pr < 3; pr++) {
	//	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	//	ind_pr[p] = 1;
	//	I->B[A_IND(I, p, i, j, k)] += 2 * I->k_viscosity_snow *
	//		(strain_rate(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], p, pr) -
	//		 strain_rate(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], p, pr)) / (2 * I->dx[pr]);
	//}
	if (p != 1) {
		ind_p[0] = ind_p[1] = ind_p[2] = 0;
		ind_p[p] = 1;
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		if (p == 0) ind_pr[2] = 1;
		else ind_pr[0] = 1;
		I->B[A_IND(I, p, i, j, k)] += 2 * I->k_viscosity_snow * (
			(velocity(I, p, i + ind_p[0], j + ind_p[1], k + ind_p[2]) - 2 * velocity(I, p, i, j, k) + velocity(I, p, i - ind_p[0], j - ind_p[1], k - ind_p[2])) / (I->dx[p] * I->dx[p]) +
			0.5 * (velocity(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - 2 * velocity(I, p, i, j, k) + velocity(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (I->dx[pr] * I->dx[pr]) +
			0.5 * (velocity(I, pr, i + ind_p[0] + ind_pr[0], j + ind_p[1] + ind_pr[1], k + ind_p[2] + ind_pr[2]) -
			velocity(I, pr, i - ind_p[0] + ind_pr[0], j - ind_p[1] + ind_pr[1], k - ind_p[2] + ind_pr[2]) -
			velocity(I, pr, i + ind_p[0] - ind_pr[0], j + ind_p[1] - ind_pr[1], k + ind_p[2] - ind_pr[2]) +
			velocity(I, pr, i - ind_p[0] - ind_pr[0], j - ind_p[1] - ind_pr[1], k - ind_p[2] - ind_pr[2])) / (I->dx[p] * I->dx[pr]));
	}
	return 0;
}

