#include "init_data.h"
#include "utils.h"
#include "x_forward_euler_second_combined_VOF.h"
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

int GRAD_pressure_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	I->B[A_IND(I, p, i, j, k)] -= (pressure(I, i + ind_p[0], j + ind_p[1], k + ind_p[2]) - pressure(I, i - ind_p[0], j - ind_p[1], k - ind_p[2])) / (2 * I->dx[p]);
	return 0;
}

int DIV_density_velocity_velocity_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / I->volume[VOLUME_IND(I, i, j, k)]) * density_on_face(I, i, j, k, s) *
				velocity_on_face(I, p, i, j, k, s) * velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_density_velocity_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / I->volume[VOLUME_IND(I, i, j, k)]) * density_on_face(I, i, j, k, s) *
				velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_linear_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] += (I->area[AREA_IND(I, i, j, k, s)] / I->volume[VOLUME_IND(I, i, j, k)]) * 2 * I->k_viscosity_snow *
				strain_rate_on_face(I, i, j, k, s, p, pr) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int VECT_barotropy_pressure_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, c1;
	c1 = 0.1;//0.001;
	A_value = 1;
	WRITE_TO_A(p, i, j, k, -1);
	A_value = - c1 * I->pressure_atmosphere / I->density_snow;
	WRITE_TO_A(3, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->pressure_atmosphere * (1 - c1);
	return 0;
}

int VECT_barotropy_density_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = 1;
	WRITE_TO_A(p, i, j, k, -1);
	A_value = - I->density_snow / I->pressure_atmosphere;
	WRITE_TO_A(4, i, j, k, -1);
}

int DIV_div_density_velocity_velocity_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, 4, i, j, k)] -= I->area[AREA_IND(I, i, j, k, s)] *
				(((density_on_face(I, i + 1, j, k, s) * velocity_on_face(I, pr, i + 1, j, k, s) * velocity_on_face(I, 0, i + 1, j, k, s) -
				   density_on_face(I, i - 1, j, k, s) * velocity_on_face(I, pr, i - 1, j, k, s) * velocity_on_face(I, 0, i - 1, j, k, s)) / (2 * I->dx[0])) +
				 ((density_on_face(I, i, j + 1, k, s) * velocity_on_face(I, pr, i, j + 1, k, s) * velocity_on_face(I, 1, i, j + 1, k, s) -
				   density_on_face(I, i, j - 1, k, s) * velocity_on_face(I, pr, i, j - 1, k, s) * velocity_on_face(I, 1, i, j - 1, k, s)) / (2 * I->dx[1])) +
				 ((density_on_face(I, i, j, k + 1, s) * velocity_on_face(I, pr, i, j, k + 1, s) * velocity_on_face(I, 2, i, j, k + 1, s) -
				   density_on_face(I, i, j, k - 1, s) * velocity_on_face(I, pr, i, j, k - 1, s) * velocity_on_face(I, 2, i, j, k - 1, s)) / (2 * I->dx[2]))) *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
}

int VECT_gravity_force_forward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	I->B[A_IND(I, p, i, j, k)] += density(I, i, j, k) * I->g[p];
	return 0;
}

