#include "init_data.h"
#include "utils.h"
#include "x_backward_euler_second_combined_VOF.h"
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

int GRAD_pressure_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int ind_p[3];
	double A_value;
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	A_value = 1 / (2 * I->dx[p]);
	WRITE_TO_A(4, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
	A_value = - A_value;
	WRITE_TO_A(4, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
	return 0;
}

int DIV_density_velocity_velocity_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * density_on_face(I, i, j, k, s) *
				velocity_on_face(I, p, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

int DIV_density_velocity_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * density_on_face(I, i, j, k, s) *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)];
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

int DIV_shear_stress_linear_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr, ind_p[3], ind_pr[3];
	double A_value;
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * I->k_viscosity_snow *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)] / (2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], s);
			A_value = - A_value;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], s);
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * I->k_viscosity_snow *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)] / (2 * I->dx[p]);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], s);
			A_value = - A_value;
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], s);
		}
	}
	return 0;
}

int VECT_barotropy_pressure_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, c1;
	c1 = 0.01;//0.001;
	A_value = 1;
	WRITE_TO_A(p, i, j, k, -1);
	A_value = - c1 * I->pressure_atmosphere / I->density_snow;
	WRITE_TO_A(3, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->pressure_atmosphere * (1 - c1);
	return 0;
}

