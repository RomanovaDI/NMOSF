#include "init_data.h"
#include "utils.h"
#include "x_backward_euler_second_separated_FDM_termogas.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

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

double divider(in *I, int p, int i, int j, int k, char object[50])
{
	int obj;
	if (strcmp(object, "density_gas_saturation_gas_porousness") == 0)
		obj = 0;
	else if (strcmp(object, "density_porousness") == 0)
		obj = 1;
	else if (strcmp(object, "porousness_density_internal_energy") == 0)
		obj = 2;
	else if (strcmp(object, "identity") == 0)
		obj = 3;
	else if (strcmp(object, "inverse_porousness_density_environment_specific_heat_environment") == 0)
		obj = 4;
	else
		printf("Error object of divider function.\n");
	double tmp;
	switch(obj) {
		case 0:
			//return 1;
			return density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->porousness;
		case 1:
			//return 1;
			return density_t(I, p - 5, i, j, k) * I->porousness;
		case 2:
			//return 1;
			tmp = 0;
			for (int pp = 0; pp < 2; pp++)
				tmp += I->porousness * density_t(I, pp, i, j, k) * saturation(I, pp, i, j, k) * I->specific_heat[pp];
			for (int pp = 0; pp < 4; pp++)
				tmp += I->porousness * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[pp + 2] * concentration(I, pp, i, j, k);
			return tmp;
		case 3:
			return 1;
		case 4:
			return (1 - I->porousness) * I->density_environment * I->specific_heat[6];
	}
}

double multiplier(in *I, int p, int pp, int pr, int i, int j, int k, char object[50])
{
	int obj;
	if (strcmp(object, "density_avarage_velocity") == 0)
		obj = 0;
	else if (strcmp(object, "density_avarage_velocity_divide_by_saturation") == 0)
		obj = 1;
	else if (strcmp(object, "density_internal_energy_avarage_velocity") == 0)
		obj = 2;
	else if (strcmp(object, "Darsi_M_coef") == 0)
		obj = 3;
	else if (strcmp(object, "thermal_conductivity") == 0)
		obj = 4;
	else if (strcmp(object, "coef_saturation") == 0)
		obj = 5;
	else if (strcmp(object, "heat_influx_vector_environment") == 0)
		obj = 6;
	else
		printf("Error object of multiplier function.\n");
	double tmp, tmp1;
	switch (obj) {
		case 0:
			return density_t(I, 2, i, j, k) * avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr]);
			//return avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, 2, i, j, k) * I->porousness);
		case 1:
			return density_t(I, p - 5, i, j, k) * avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr]);
			return avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, p - 5, i, j, k) * I->porousness);
			//return -density_t(I, p - 5, i, j, k) * I->permeability * relative_permeability_derivative_with_recpect_to_saturation(I, p - 5, pp - 5, i, j, k) * grad_pressure(I, pr, i, j, k) /
			//		(2 * I->dx[pr] * viscosity(I, p - 5, i, j, k));
		case 2:
			//return avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr] * I->porousness);
			tmp = 0;
			for (pp = 0; pp < 2; pp++)
				tmp += density_t(I, pp, i, j, k) * I->specific_heat[pp] * I->adiabatic_exponent[pp] * avarage_velocity(I, pp, pr, i, j, k) / (2 * I->dx[pr]);
			for (pp = 0; pp < 4; pp++)
				tmp += density_t(I, 2, i, j, k) * I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] * concentration(I, pp, i, j, k) *
					avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr]);
			return tmp;
			/*tmp1 = 0;
			for (pp = 0; pp < 2; pp++)
				tmp1 += I->porousness * density_t(I, pp, i, j, k) * saturation(I, pp, i, j, k) * I->specific_heat[pp];
			for (pp = 0; pp < 4; pp++)
				tmp1 += I->porousness * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[pp + 2] * concentration(I, pp, i, j, k);
			return tmp / tmp1;*/
		case 3:
			return Darsi_M_coef(I, i, j, k) / (2 * I->dx[pr] * I->dx[pr]); 
		case 4:
			return I->porousness * thermal_conductivity(I, i, j, k) / (2 * I->dx[pr] * I->dx[pr]);
		case 5:
			return density_t(I, p - 5, i, j, k) * I->permeability * relative_permeability(I, p - 5, i, j, k) * coef_grad_saturation(I, p - 5, pp - 5, i, j, k) *
					capillary_pressure_derivative_by_saturation(I, 0, i, j, k) / viscosity(I, p - 5, i, j, k);
		case 6:
			return (1 - I->porousness) * I->thermal_conductivity_coef[3] / (2 * I->dx[pr] * I->dx[pr]);

	}
}

double operand(in * I, int p, int i, int j, int k, char object[50])
{
	int obj;
	if (strcmp(object, "concentration") == 0)
		obj = 0;
	else if (strcmp(object, "saturation") == 0)
		obj = 1;
	else if (strcmp(object, "temperature_flow") == 0)
		obj = 2;
	else if (strcmp(object, "pressure") == 0)
		obj = 3;
	else if (strcmp(object, "temperature_environment") == 0)
		obj = 4;
	else if (strcmp(object, "temperature_flow_minus_initial") == 0)
		obj = 5;
	else if (strcmp(object, "temperature_environment_minus_initial") == 0)
		obj = 6;
	else
		printf("Error object of operand function.\n");
	switch (obj){
		case 0:
			return concentration(I, p, i, j, k);
		case 1:
			return saturation(I, p - 5, i, j, k);
		case 2:
			return temperature_flow(I, i, j, k);
		case 3:
			return pressure(I, i, j, k);
		case 4:
			return temperature_environment(I, i, j, k);
		case 5:
			return temperature_flow(I, i, j, k) - I->initial_temperature;
		case 6:
			return temperature_environment(I, i, j, k) - I->initial_temperature;
	}
}
/*
double lambda2D_laplas(in *i, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50])
{
	double lambda = 0;
	for (int pr = 0; pr < 2; pr ++) {
		if ()
	}
}
*/
double lambda2D(in *I, int p, int pp, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_divider[50])
{
	double lambda = 0, omega = 0, sign = 0;
	for (int pr = 0; pr < 2; pr ++) {
		if ((pr == 0) && (i == ii) && (j == jj) && !(boundary_cell(I, i + 1, j, k)) && !(boundary_cell(I, i - 1, j, k))) {
			for (int jjj = j - 1; jjj <= j + 1; jjj++) {
				if (jjj == j)
					omega = 0.75;
				else
					omega = 0.125;
				lambda += omega * (multiplier(I, p, pp, pr, ii + 1, jjj, kk, object_multiplier) - multiplier(I, p, pp, pr, ii - 1, jjj, kk, object_multiplier));
			}
		} else if ((pr == 1) && (i == ii) && (j == jj) && !(boundary_cell(I, i, j + 1, k)) && !(boundary_cell(I, i, j - 1, k))) {
			for (int iii = i - 1; iii <= i + 1; iii++) {
				if (iii == i)
					omega = 0.75;
				else
					omega = 0.125;
				lambda += omega * (multiplier(I, p, pp, pr, iii, jj + 1, kk, object_multiplier) - multiplier(I, p, pp, pr, iii, jj - 1, kk, object_multiplier));
			}
		} else {
			if ((((pr == 0) && (j == jj)) || ((pr == 1) && (i == ii))))
				omega = 0.75;
			else
				omega = 0.125;
			if ((pr == 0) && (ii > i))
				sign = -1;
			else if ((pr == 0) && (ii < i))
				sign = 1;
			else if ((pr == 0) && (i == ii))
				sign = 0;
			else if ((pr == 1) && (jj > j))
				sign = -1;
			else if ((pr == 1) && (jj < j))
				sign = 1;
			else if ((pr == 1) && (j == jj))
				sign = 0;
			lambda += sign * omega * multiplier(I, p, pp, pr, ii, jj, kk, object_multiplier);
		}
	}
	lambda /= divider(I, p, i, j, k, object_divider);
	return lambda;
}

double delta(in *I, int p, int pp, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_divider[50])
{
	if ((i == ii) && (j == jj) && (k == kk))
		return 0;
	else
		return max3(0, - lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider), - lambda2D(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider));
}

double alpha(in *I, int p, int pp, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_operand[50], char object_divider[50])
{
	if ((i == ii) && (j == jj) && (k == kk))
		return 0;
	double Q1 = 0, Q2 = 0, P1 = 0, P2 = 0;
	int kkk = k;
	for (int iii = i - 1; iii < i + 2; iii++) {
		for (int jjj = j - 1; jjj < j + 2; jjj++) {
			if (!((iii == 0) && (jjj == 0))) {
				Q1 += max2(0, lambda2D(I, p, pp, i, j, k, iii, jjj, kkk, object_multiplier, object_divider)) * max2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, i, j, k, object_operand), 0);
				Q2 += max2(0, lambda2D(I, p, pp, i, j, k, iii, jjj, kkk, object_multiplier, object_divider)) * min2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, i, j, k, object_operand), 0);
				P1 += min2(0, lambda2D(I, p, pp, i, j, k, iii, jjj, kkk, object_multiplier, object_divider)) * min2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, i, j, k, object_operand), 0);
				P2 += min2(0, lambda2D(I, p, pp, i, j, k, iii, jjj, kkk, object_multiplier, object_divider)) * max2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, i, j, k, object_operand), 0);
			}
		}
	}
	if (P1 == 0)
		Q1 = 0;
	else
		Q1 /= P1;
	if (P2 == 0)
		Q2 = 0;
	else
		Q2 /= P2;
	double R1 = max2(0, min2(1, Q1));
	double R2 = max2(0, min2(1, Q2));
	if (lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) <= lambda2D(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider)) {
		if (operand(I, pp, i, j, k, object_operand) >= operand(I, pp, ii, jj, kk, object_operand))
			return min2(R1 * delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider), lambda2D(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider) + delta(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider));
		else
			return min2(R2 * delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider), lambda2D(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider) + delta(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider));
	} else {
		Q1 = Q2 = P1 = P2 = 0;
		int kkk = kk;
		for (int iii = ii - 1; iii < ii + 2; iii++) {
			for (int jjj = jj - 1; jjj < jj + 2; jjj++) {
				if (!((iii == 0) && (jjj == 0))) {
					Q1 += max2(0, lambda2D(I, p, pp, ii, jj, kk, iii, jjj, kkk, object_multiplier, object_divider)) * max2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, ii, jj, kk, object_operand), 0);
					Q2 += max2(0, lambda2D(I, p, pp, ii, jj, kk, iii, jjj, kkk, object_multiplier, object_divider)) * min2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, ii, jj, kk, object_operand), 0);
					P1 += min2(0, lambda2D(I, p, pp, ii, jj, kk, iii, jjj, kkk, object_multiplier, object_divider)) * min2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, ii, jj, kk, object_operand), 0);
					P2 += min2(0, lambda2D(I, p, pp, ii, jj, kk, iii, jjj, kkk, object_multiplier, object_divider)) * max2(operand(I, pp, iii, jjj, kkk, object_operand) - operand(I, pp, ii, jj, kk, object_operand), 0);
				}
			}
		}
		if (P1 == 0)
			Q1 = 0;
		else
			Q1 /= P1;
		if (P2 == 0)
			Q2 = 0;
		else
			Q2 /= P2;
		double R1 = max2(0, min2(1, Q1));
		double R2 = max2(0, min2(1, Q2));
		if (operand(I, pp, ii, jj, kk, object_operand) >= operand(I, pp, i, j, k, object_operand))
			return min2(R1 * delta(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider), lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) + delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider));
		else
			return min2(R2 * delta(I, p, pp, ii, jj, kk, i, j, k, object_multiplier, object_divider), lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) + delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider));
	}
}

double lambda_res(in *I, int p, int pp, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_operand[50], char object_divider[50], int antidiff)
{
	if ((i == ii) && (j == jj) && (k == kk))
		return lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider);
	else if (antidiff == 2)
		return lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) + delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) - alpha(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_operand, object_divider);
	else if (antidiff == 1)
		return lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider) + delta(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider);
	else if (antidiff == 0)
		return lambda2D(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_divider);
	else
		printf("Error value of antidiff.\n");
	return 0;
}

int div_func(in *I, int p, int pp, int i, int j, int k, char object_multiplier[50], char object_operand[50], char object_divider[50], int antidiff, double tetta)
{
	if (k != 0)
		return 0;
	for (int ii = i - 1; ii < i + 2; ii++) {
		for (int jj = j - 1; jj < j + 2; jj++) {
			double A_value, tmp;
			int kk = k;
			if (!boundary_cell(I, ii, jj, kk)) {
			//if (1) {
				tmp = lambda_res(I, p, pp, i, j, k, ii, jj, kk, object_multiplier, object_operand, object_divider, antidiff);
				A_value = tetta * tmp;
				WRITE_TO_A(pp, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= (1 - tetta) * tmp * operand(I, pp, i, j, k, object_operand);
				if (!((i == ii) && (j == jj) && (k == kk))) {
					A_value *= -1;
					WRITE_TO_A(pp, ii, jj, kk, -1);
					I->B[A_IND(I, p, i, j, k)] += (1 - tetta) * tmp * operand(I, pp, ii, jj, kk, object_operand);
				}
			}
		}
	}
	return 0;
}

int lapl_func(in *I, int p, int pp, int i, int j, int k, char object_multiplier[50], char object_operand[50], char object_divider[50])
{
	for (int pr = 0; pr < 3; pr++) {
		double A_value;
		int ind_pr[3];
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))) {
			A_value = (multiplier(I, p, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], object_multiplier) + multiplier(I, p, pp, pr, i, j, k, object_multiplier)) / divider(I, p, i, j, k, object_divider);
			WRITE_TO_A(pp, i, j, k, -1);
			A_value = - A_value;
			WRITE_TO_A(pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		}
		if (!(boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			A_value = (multiplier(I, p, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], object_multiplier) + multiplier(I, p, pp, pr, i, j, k, object_multiplier)) / divider(I, p, i, j, k, object_divider);
			WRITE_TO_A(pp, i, j, k, -1);
			A_value = - A_value;
			WRITE_TO_A(pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

int DIV_concentration_density_average_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 0) && (p != 1) && (p != 2) && (p != 3)) {
		printf("Incorrect index of DIV_concentration_density_average_velocity\n");
		return 1;
	}
	if (div_func(I, p, p, i, j, k, "density_avarage_velocity", "concentration", "density_gas_saturation_gas_porousness", 1, 0.5)) return 1;
	return 0;
}

int SCAL_mass_inflow_rate_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, tetta = 1;
	/*if ((p == 6) || (p == 7)) {
		A_value = - mass_inflow_rate_func(I, p, i, j, k) / saturation(I, p - 5, i, j, k);
		WRITE_TO_A(p, i, j, k, -1);
	} else if (p == 1) {
		A_value = - mass_inflow_rate_func(I, p, i, j, k) / concentration(I, p, i, j, k);
		WRITE_TO_A(p, i, j, k, -1);
	} else {
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k);
	}*/
/*
	if ((p == 0) || (p == 1) || (p == 2) || (p == 3))
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k) / divider(I, p, i, j, k, "density_gas_saturation_gas_porousness");
	else if ((p == 5) || (p == 6) || (p == 7))
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k) / divider(I, p, i, j, k, "density_porousness");
	else
		printf("Error SCAL_mass_inflow_rate_backward_euler_second_separated_FDM_termogas index.\n");
*/
	if ((p == 0) || (p == 1) || (p == 2) || (p == 3)) {
		I->B[A_IND(I, p, i, j, k)] += (1 - tetta) * mass_inflow_rate_func(I, p, i, j, k) / divider(I, p, i, j, k, "density_gas_saturation_gas_porousness");
		A_value = - tetta * mass_inflow_rate_func(I, p, i, j, k) / (divider(I, p, i, j, k, "density_gas_saturation_gas_porousness") * I->B_prev[B_IND(I, p, i, j, k)]);
		WRITE_TO_A(p, i, j, k, -1);
	} else if ((p == 5) || (p == 6) || (p == 7)) {
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k) / divider(I, p, i, j, k, "density_porousness");
		A_value = - tetta * mass_inflow_rate_func(I, p, i, j, k) / (divider(I, p, i, j, k, "density_porousness") * I->B_prev[B_IND(I, p, i, j, k)]);
		WRITE_TO_A(p, i, j, k, -1);
	} else
		printf("Error SCAL_mass_inflow_rate_backward_euler_second_separated_FDM_termogas index.\n");
	return 0;
}

int LAPL_coef_pressure_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (p != 4) {
		printf("Error LAPL_coef_pressure_backward_euler_second_separated_FDM_termogas index.\n");
	}
	if (lapl_func(I, p, p, i, j, k, "Darsi_M_coef", "pressure", "identity")) return 1;
	return 0;
}

int DIV_density_average_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 5) && (p != 6) && (p != 7)) {
		printf("Incorrect index of DIV_density_average_velocity\n");
		return 1;
	}
	if (div_func(I, p, p, i, j, k, "density_avarage_velocity_divide_by_saturation", "saturation", "density_porousness", 1, 0.5)) return 1;
	/*for (int pp = 5; pp < 8; pp++) {
		if (div_func(I, p, pp, i, j, k, "density_avarage_velocity_divide_by_saturation", "saturation", "density_porousness", 1, 0.5)) return 1;
		for (int pr = 0; pr < 3; pr++) {
			int ind_pr[3];
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			I->B[A_IND(I, p, i, j, k)] += (multiplier(I, p, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], "density_avarage_velocity_divide_by_saturation") -
				multiplier(I, p, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], "density_avarage_velocity_divide_by_saturation")) * saturation(I, pp - 5, i, j, k);
		}
	}*/
	//if (LAPL(p, i, j, k, coef_saturation, backward_euler, second, separated, FDM, termogas)) return 1;
	return 0;
}

int LAPL_coef_saturation_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 5) && (p != 6) && (p != 7)) {
		printf("Error LAPL_coef_saturation_backward_euler_second_separated_FDM_termogas index.\n");
	}
	if (lapl_func(I, p, 5, i, j, k, "coef_saturation", "saturation", "density_porousness")) return 1;
	if (lapl_func(I, p, 7, i, j, k, "coef_saturation", "saturation", "density_porousness")) return 1;
	return 0;
}

int DIV_density_internal_energy_avarage_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (p != 8) {
		printf("Incorrect index of DIV_density_internal_energy_average_velocity\n");
		return 1;
	}
	if (div_func(I, p, p, i, j, k, "density_internal_energy_avarage_velocity", "temperature_flow_minus_initial", "porousness_density_internal_energy", 1, 0.5)) return 1;
	//if (div_func(I, p, p, i, j, k, "density_internal_energy_avarage_velocity", "temperature_flow", "porousness_density_internal_energy", 1, 0.5)) return 1;
	//if (div_func(I, p, p, i, j, k, "density_internal_energy_avarage_velocity", "temperature_flow", "identity", 1, 0.5)) return 1;
	return 0;
}

int LAPL_heat_influx_vector_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (lapl_func(I, p, p, i, j, k, "thermal_conductivity", "temperature_flow", "porousness_density_internal_energy")) return 1;
	return 0;
}

int SCAL_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = I->heat_transfer_coef / divider(I, p, i, j, k, "porousness_density_internal_energy");
	WRITE_TO_A(p, i, j, k, -1);
	//A_value *= -1;
	//WRITE_TO_A(9, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->heat_transfer_coef * (temperature_environment(I, i, j, k) - I->initial_temperature) / divider(I, p, i, j, k, "porousness_density_internal_energy");
	return 0;
}

int SCAL_chemical_reaction_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
/*
	for (int pp = 0; pp < 4; pp++) {
		I->B[A_IND(I, p, i, j, k)] -= mass_inflow_rate_func(I, pp, i, j, k) * I->initial_enthalpy[pp + 2] / divider(I, p, i, j, k, "porousness_density_internal_energy");
	}
	for (int pp = 0; pp < 2; pp++) {
		I->B[A_IND(I, p, i, j, k)] -= mass_inflow_rate_func(I, pp + 5, i, j, k) * I->initial_enthalpy[pp] / divider(I, p, i, j, k, "porousness_density_internal_energy");
	}
*/
	double A_value, tetta = 0.5;
	for (int pp = 0; pp < 4; pp++) {
		A_value += tetta * mass_inflow_rate_func(I, pp, i, j, k) * I->initial_enthalpy[pp + 2] /
			(divider(I, p, i, j, k, "porousness_density_internal_energy") * temperature_flow(I, i, j, k));
		I->B[A_IND(I, p, i, j, k)] -= (1 - tetta) * mass_inflow_rate_func(I, pp, i, j, k) * I->initial_enthalpy[pp + 2] / divider(I, p, i, j, k, "porousness_density_internal_energy");
	}
	for (int pp = 0; pp < 2; pp++) {
		A_value += tetta * mass_inflow_rate_func(I, pp + 5, i, j, k) * I->initial_enthalpy[pp] /
			(divider(I, p, i, j, k, "porousness_density_internal_energy") * temperature_flow(I, i, j, k));
		I->B[A_IND(I, p, i, j, k)] -= (1 - tetta) * mass_inflow_rate_func(I, pp + 5, i, j, k) * I->initial_enthalpy[pp] / divider(I, p, i, j, k, "porousness_density_internal_energy");
	}
	WRITE_TO_A(p, i, j, k, -1);
	return 0;
}

int LAPL_heat_influx_vector_environment_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (lapl_func(I, p, p, i, j, k, "heat_influx_vector_environment", "temperature_environment_minus_initial", "inverse_porousness_density_environment_specific_heat_environment")) return 1;
	return 0;
}

int SCAL_minus_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = I->heat_transfer_coef / divider(I, p, i, j, k, "inverse_porousness_density_environment_specific_heat_environment");
	WRITE_TO_A(p, i, j, k, -1);
	//A_value *= -1;
	//WRITE_TO_A(8, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->heat_transfer_coef * (temperature_flow(I, i, j, k) - I->initial_temperature) / divider(I, p, i, j, k, "inverse_porousness_density_environment_specific_heat_environment");
	return 0;
}

