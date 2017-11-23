#include "init_data.h"
#include "utils.h"
#include "x_crank_nikolson_second_combined_VOF_avalanche.h"
#include "mpi.h"
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

int DIV_density_velocity_velocity_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * density_on_face(I, i, j, k, s) *
				velocity_on_face(I, p, i, j, k, s) * velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_density_velocity_velocity_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
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

int DIV_density_velocity_velocity_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int GRAD_pressure_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	I->B[A_IND(I, p, i, j, k)] -= (pressure(I, i + ind_p[0], j + ind_p[1], k + ind_p[2]) - pressure(I, i - ind_p[0], j - ind_p[1], k - ind_p[2])) / (4 * I->dx[p]);
	return 0;
}

int GRAD_pressure_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	A_value = 1 / (4 * I->dx[p]);
	WRITE_TO_A(4, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
	A_value = - A_value;
	WRITE_TO_A(4, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
	return 0;
}

int GRAD_pressure_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (GRAD(p, i, j, k, pressure, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int VECT_gravity_force_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	I->B[A_IND(I, p, i, j, k)] += density(I, i, j, k) * I->g[p] / 2;
	return 0;
}

int VECT_gravity_force_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	double A_value;
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	A_value = - (I->density_snow - I->density_air) * I->g[p] / 2;
	WRITE_TO_A(3, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->density_air * I->g[p] / 2;
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_shear_stress_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] += (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * effective_viscosity_on_face(I, i, j, k, s) *
				strain_rate_on_face(I, i, j, k, s, p, pr) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr, ind_p[3], ind_pr[3];
	double d_p, d_pr, A_value;
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	d_p = I->dx[p];
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			d_pr = I->dx[pr];
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * effective_viscosity_on_face(I, i, j, k, s) *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)] / (4 * d_pr);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], s);
			A_value = - A_value;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], s);
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * effective_viscosity_on_face(I, i, j, k, s) *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)] / (4 * d_p);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], s);
			A_value = - A_value;
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], s);
		}
	}
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_grad_pressure_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		I->B[A_IND(I, 4, i, j, k)] -=  (I->area[AREA_IND(I, i, j, k, s)] / 2) * (
				(pressure_on_face(I, i + 1, j, k, s) - pressure_on_face(I, i - 1, j, k, s)) * I->normal[NORMAL_IND(I, 0, i, j, k, s)] / (2 * I->dx[0]) +
				(pressure_on_face(I, i, j + 1, k, s) - pressure_on_face(I, i, j - 1, k, s)) * I->normal[NORMAL_IND(I, 1, i, j, k, s)] / (2 * I->dx[1]) +
				(pressure_on_face(I, i, j, k + 1, s) - pressure_on_face(I, i, j, k - 1, s)) * I->normal[NORMAL_IND(I, 2, i, j, k, s)] / (2 * I->dx[2]));
	}
	return 0;
}

int DIV_grad_pressure_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 4.) * I->normal[NORMAL_IND(I, 0, i, j, k, s)] / (2. * I->dx[0]);
		WRITE_TO_A(4, i + 1, j, k, -1);
		WRITE_TO_A(4, i + 1, j, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i - 1, j, k, -1);
		WRITE_TO_A(4, i - 1, j, k, s);
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 4.) * I->normal[NORMAL_IND(I, 1, i, j, k, s)] / (2. * I->dx[1]);
		WRITE_TO_A(4, i, j + 1, k, -1);
		WRITE_TO_A(4, i, j + 1, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j - 1, k, -1);
		WRITE_TO_A(4, i, j - 1, k, s);
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 4.) * I->normal[NORMAL_IND(I, 2, i, j, k, s)] / (2. * I->dx[2]);
		WRITE_TO_A(4, i, j, k + 1, -1);
		WRITE_TO_A(4, i, j, k + 1, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j, k - 1, -1);
		WRITE_TO_A(4, i, j, k - 1, s);
	}
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, 4, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / 2) *
				(((density_on_face(I, i + 1, j, k, s) * velocity_on_face(I, pr, i + 1, j, k, s) * velocity_on_face(I, 0, i + 1, j, k, s) -
				   density_on_face(I, i - 1, j, k, s) * velocity_on_face(I, pr, i - 1, j, k, s) * velocity_on_face(I, 0, i - 1, j, k, s)) / (2 * I->dx[0])) +
				 ((density_on_face(I, i, j + 1, k, s) * velocity_on_face(I, pr, i, j + 1, k, s) * velocity_on_face(I, 1, i, j + 1, k, s) -
				   density_on_face(I, i, j - 1, k, s) * velocity_on_face(I, pr, i, j - 1, k, s) * velocity_on_face(I, 1, i, j - 1, k, s)) / (2 * I->dx[1])) +
				 ((density_on_face(I, i, j, k + 1, s) * velocity_on_face(I, pr, i, j, k + 1, s) * velocity_on_face(I, 2, i, j, k + 1, s) -
				   density_on_face(I, i, j, k - 1, s) * velocity_on_face(I, pr, i, j, k - 1, s) * velocity_on_face(I, 2, i, j, k - 1, s)) / (2 * I->dx[2]))) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr, pp;
	double A_value, d_pr, ind_pr[3];
	for (s = 0; s < 6; s++) {
		for (pp = 0; pp < 3; pp++) {
			for (pr = 0; pr < 3; pr++) {
				ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
				ind_pr[pr] = 1;
				d_pr = I->dx[pr];
				A_value = (I->area[AREA_IND(I, i, j, k, s)] / 4) * density_on_face(I, i, j, k, s) * velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pp, i, j, k, s)] / (2 * d_pr);
				WRITE_TO_A(pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
				WRITE_TO_A(pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], s);
				A_value = - A_value;
				WRITE_TO_A(pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
				WRITE_TO_A(pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], s);
			}
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_velocity_cont_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, s;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] -= 0.5 * velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_velocity_cont_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, s;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = 0.5 * 0.5 * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

int DIV_velocity_cont_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, velocity_cont, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, velocity_cont, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_snow_volume_fraction_velocity_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
	I->B[A_IND(I, 3, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * phase_fraction_on_face(I, i, j, k, s) *
		(velocity_on_face(I, 0, i, j, k, s) * I->normal[NORMAL_IND(I, 0, i, j, k, s)] +
		velocity_on_face(I, 1, i, j, k, s) * I->normal[NORMAL_IND(I, 1, i, j, k, s)] +
		velocity_on_face(I, 2, i, j, k, s) * I->normal[NORMAL_IND(I, 2, i, j, k, s)]);
	}
	return 0;
}

int DIV_snow_volume_fraction_velocity_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
	A_value = (I->area[AREA_IND(I, i ,j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) *
		(velocity_on_face(I, 0, i, j, k, s) * I->normal[NORMAL_IND(I, 0, i, j, k, s)] +
		velocity_on_face(I, 1, i, j, k, s) * I->normal[NORMAL_IND(I, 1, i, j, k, s)] +
		velocity_on_face(I, 2, i, j, k, s) * I->normal[NORMAL_IND(I, 2, i, j, k, s)]);
	WRITE_TO_A(3, i, j, k, -1);
	WRITE_TO_A(3, i, j, k, s);
	}
	return 0;
}

int DIV_snow_volume_fraction_velocity_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, snow_volume_fraction_velocity, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, snow_volume_fraction_velocity, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_density_velocity_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] -= (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * density_on_face(I, i, j, k, s) *
				velocity_on_face(I, pr, i, j, k, s) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_density_velocity_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = (I->area[AREA_IND(I, i, j, k, s)] / (4 * I->volume[VOLUME_IND(I, i, j, k)])) * density_on_face(I, i, j, k, s) *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)];
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

int DIV_density_velocity_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, density_velocity, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int DIV_shear_stress_linear_half_forward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			I->B[A_IND(I, p, i, j, k)] += (I->area[AREA_IND(I, i, j, k, s)] / (2 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * I->k_viscosity_snow *
				strain_rate_on_face(I, i, j, k, s, p, pr) * I->normal[NORMAL_IND(I, pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_linear_half_backward_euler_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
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
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (8 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * I->k_viscosity_snow *
				I->normal[NORMAL_IND(I, pr, i, j, k, s)] / (2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], s);
			A_value = - A_value;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], s);
			A_value = - (I->area[AREA_IND(I, i, j, k, s)] / (8 * I->volume[VOLUME_IND(I, i, j, k)])) * 2 * I->k_viscosity_snow *
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

int DIV_shear_stress_linear_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress_linear, half_forward_euler, second, combined, VOF, avalanche)) return 1;
	if (DIV(p, i, j, k, shear_stress_linear, half_backward_euler, second, combined, VOF, avalanche)) return 1;
	return 0;
}

int VECT_barotropy_pressure_crank_nikolson_second_combined_VOF_avalanche(in *I, int p, int i, int j, int k)
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

