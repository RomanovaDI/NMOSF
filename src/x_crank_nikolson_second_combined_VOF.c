#include "utils.h"
#include "x_crank_nikolson_second_combined.h"

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DIV_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int GRAD_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int VECT_gravity_force_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_shear_stress_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_grad_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_density_snow_volume_fraction_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_velocity_cont_crank_nikolson_second_combined_VOF(int p, int i, int j, int k);
int DIV_velocity_cont_half_forward_euler_second_combined_VOF(int p, int i, int j, int k);
int DIV_velocity_cont_half_backward_euler_second_combined_VOF(int p, int i, int j, int k);

int DIV_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, density_velocity_velocity, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] -= (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * density_on_face(i, j, k, s) *
				velocity_on_face(p, i, j, k, s) * velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * density_on_face(i, j, k, s) *
			(velocity_on_face(0, i, j, k, s) * normal[NORMAL_IND(0, i, j, k, s)] +
			 velocity_on_face(1, i, j, k, s) * normal[NORMAL_IND(1, i, j, k, s)] +
			 velocity_on_face(2, i, j, k, s) * normal[NORMAL_IND(2, i, j, k, s)]);
		WRITE_TO_A(p, i, j, k, -1);
		WRITE_TO_A(p, i, j, k, s);
	}
	return 0;
}

int GRAD_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (GRAD(p, i, j, k, pressure, half_forward_euler, second, combined, VOF)) return 1;
	if (GRAD(p, i, j, k, pressure, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int GRAD_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	B[A_IND(p, i, j, k)] -= (pressure(i + ind_p[0], j + ind_p[1], k + ind_p[2]) - pressure(i - ind_p[0], j - ind_p[1], k - ind_p[2])) / dx[p];
	return 0;
}

int GRAD_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	int ind_p[3];
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	A_value = 1 / (4 * dx[p]);
	WRITE_TO_A(4, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
	A_value = - A_value;
	WRITE_TO_A(4, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
	return 0;
}

int VECT_gravity_force_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_forward_euler, second, combined, VOF)) return 1;
	if (VECT(p, i, j, k, gravity_force, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int VECT_gravity_force_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	B[A_IND(p, i, j, k)] += density(i, j, k) * g[p] / 2;
	return 0;
}

int VECT_gravity_force_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	double A_value;
	if (check_for_corrupt_cell(i, j, k)) return 1;
	A_value = - (density_snow - density_air) * g[p] / 2;
	WRITE_TO_A(3, i, j, k, -1);
	B[A_IND(p, i, j, k)] += density_air * g[p] / 2;
	return 0;
}

int DIV_shear_stress_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, shear_stress, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_shear_stress_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] += (area[AREA_IND(i, j, k, s)] / (2 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				strain_rate_on_face(i, j, k, s, p, pr) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_shear_stress_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, ind_p[3], ind_pr[3];
	double d_p, d_pr, A_value;
	ind_p[0] = ind_p[1] = ind_p[2] = 0;
	ind_p[p] = 1;
	d_p = dx[p];
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			d_pr = dx[pr];
			A_value = - (area[AREA_IND(i, j, k, s)] / (4 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				normal[NORMAL_IND(pr, i, j, k, s)] / (4 * d_pr);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], s);
			A_value = - A_value;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], s);
			A_value = - (area[AREA_IND(i, j, k, s)] / (4 * volume[VOLUME_IND(i, j, k)])) * 2 * effective_viscosity_on_face(i, j, k, s) *
				normal[NORMAL_IND(pr, i, j, k, s)] / (4 * d_p);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], -1);
			WRITE_TO_A(pr, i + ind_p[0], j + ind_p[1], k + ind_p[2], s);
			A_value = - A_value;
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], -1);
			WRITE_TO_A(pr, i - ind_p[0], j - ind_p[1], k - ind_p[2], s);
		}
	}
	return 0;
}

int DIV_grad_pressure_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, grad_pressure, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_grad_pressure_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	for (s = 0; s < 6; s++) {
		B[A_IND(4, i, j, k)] -=  (area[AREA_IND(i, j, k, s)] / 2) * (
				(pressure_on_face(i + 1, j, k, s) - pressure_on_face(i - 1, j, k, s)) * normal[NORMAL_IND(0, i, j, k, s)] / (2 * dx[0]) +
				(pressure_on_face(i, j + 1, k, s) - pressure_on_face(i, j - 1, k, s)) * normal[NORMAL_IND(1, i, j, k, s)] / (2 * dx[1]) +
				(pressure_on_face(i, j, k + 1, s) - pressure_on_face(i, j, k - 1, s)) * normal[NORMAL_IND(2, i, j, k, s)] / (2 * dx[2]));
	}
	return 0;
}

int DIV_grad_pressure_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (area[AREA_IND(i, j, k, s)] / 4.) * normal[NORMAL_IND(0, i, j, k, s)] / (2. * dx[0]);
		WRITE_TO_A(4, i + 1, j, k, -1);
		WRITE_TO_A(4, i + 1, j, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i - 1, j, k, -1);
		WRITE_TO_A(4, i - 1, j, k, s);
		A_value = (area[AREA_IND(i, j, k, s)] / 4.) * normal[NORMAL_IND(1, i, j, k, s)] / (2. * dx[1]);
		WRITE_TO_A(4, i, j + 1, k, -1);
		WRITE_TO_A(4, i, j + 1, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j - 1, k, -1);
		WRITE_TO_A(4, i, j - 1, k, s);
		A_value = (area[AREA_IND(i, j, k, s)] / 4.) * normal[NORMAL_IND(2, i, j, k, s)] / (2. * dx[2]);
		WRITE_TO_A(4, i, j, k + 1, -1);
		WRITE_TO_A(4, i, j, k + 1, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j, k - 1, -1);
		WRITE_TO_A(4, i, j, k - 1, s);
	}
	return 0;
}

int DIV_div_density_velocity_velocity_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, div_density_velocity_velocity, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_div_density_velocity_velocity_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(4, i, j, k)] -= (area[AREA_IND(i, j, k, s)] / 2) *
				(((density_on_face(i + 1, j, k, s) * velocity_on_face(pr, i + 1, j, k, s) * velocity_on_face(0, i + 1, j, k, s) -
				   density_on_face(i - 1, j, k, s) * velocity_on_face(pr, i - 1, j, k, s) * velocity_on_face(0, i - 1, j, k, s)) / (2 * dx[0])) +
				 ((density_on_face(i, j + 1, k, s) * velocity_on_face(pr, i, j + 1, k, s) * velocity_on_face(1, i, j + 1, k, s) -
				   density_on_face(i, j - 1, k, s) * velocity_on_face(pr, i, j - 1, k, s) * velocity_on_face(1, i, j - 1, k, s)) / (2 * dx[1])) +
				 ((density_on_face(i, j, k + 1, s) * velocity_on_face(pr, i, j, k + 1, s) * velocity_on_face(2, i, j, k + 1, s) -
				   density_on_face(i, j, k - 1, s) * velocity_on_face(pr, i, j, k - 1, s) * velocity_on_face(2, i, j, k - 1, s)) / (2 * dx[2]))) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_div_density_velocity_velocity_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int s, pr, pp;
	double A_value, d_pr, ind_pr[3];
	for (s = 0; s < 6; s++) {
		for (pp = 0; pp < 3; pp++) {
			for (pr = 0; pr < 3; pr++) {
				ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
				ind_pr[pr] = 1;
				d_pr = dx[pr];
				A_value = (area[AREA_IND(i, j, k, s)] / 4) * density_on_face(i, j, k, s) * velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pp, i, j, k, s)] / (2 * d_pr);
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

int DIV_velocity_cont_crank_nikolson_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	if (DIV(p, i, j, k, velocity_cont, half_forward_euler, second, combined, VOF)) return 1;
	if (DIV(p, i, j, k, velocity_cont, half_backward_euler, second, combined, VOF)) return 1;
	return 0;
}

int DIV_velocity_cont_half_forward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pr, s;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			B[A_IND(p, i, j, k)] -= 0.5 * velocity_on_face(pr, i, j, k, s) * normal[NORMAL_IND(pr, i, j, k, s)];
		}
	}
	return 0;
}

int DIV_velocity_cont_half_backward_euler_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	int pr, s;
	double A_value;
	for (s = 0; s < 6; s++) {
		for (pr = 0; pr < 3; pr++) {
			A_value = 0.5 * 0.5 * normal[NORMAL_IND(pr, i, j, k, s)];
			WRITE_TO_A(pr, i, j, k, -1);
			WRITE_TO_A(pr, i, j, k, s);
		}
	}
	return 0;
}

