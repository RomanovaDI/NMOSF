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

double multiplier(in *I, int p, int pr, int i, int j, int k, char object[50])
{
	int obj;
	if (strcmp(object, "density_avarage_velocity") == 0)
		obj = 0;
	else if (strcmp(object, "density_avarage_velocity_divide_by_saturation") == 0)
		obj = 1;
	else if (strcmp(object, "density_internal_energy_avarage_velocity") == 0)
		obj = 2;
	double tmp = 0;
	switch (obj) {
		case 0:
			return density_t(I, 2, i, j, k) * avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr]);
			//return avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, 2, i, j, k) * I->porousness);
			//return density_t(I, 2, i, j, k) * avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr] * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->porousness);
		case 1:
			//return density_t(I, p - 5, i, j, k) * avarage_velocity_derivative_with_respect_to_saturation(I, p - 5, pr, i, j, k) / (2 * I->dx[pr]);
			return density_t(I, p - 5, i, j, k) * avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, p - 5, i, j, k));
			//return avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr] * I->porousness);
			//return avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, p - 5, i, j, k) * I->porousness);
			//return density_t(I, p - 5, i, j, k) * avarage_velocity(I, p - 5, pr, i, j, k) / (2 * I->dx[pr] * saturation(I, p - 5, i, j, k) * density_t(I, p - 5, i, j, k) * I->porousness);
		case 2:
			tmp = 0;
			for (int pp = 0; pp < 1; pp++)
				tmp += density_t(I, pp, i, j, k) *
					I->specific_heat[pp] * I->adiabatic_exponent[pp] *
					avarage_velocity(I, pp, pr, i, j, k) / (2 * I->dx[pr]);
			for (int pp = 0; pp < 0; pp++)
				tmp += density_t(I, 2, i, j, k) *
					I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] *
					concentration(I, pp, i, j, k) *
					avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr]);
			return tmp;
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
	switch (obj){
		case 0:
			return concentration(I, p, i, j, k);
		case 1:
			//return 1;
			return saturation(I, p - 5, i, j, k);
		case 2:
			return temperature_flow(I, i, j, k);
	}
}

double lambda2D(in *I, int p, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50])
{
	double lambda = 0, omega = 0, sign = 0;
	for (int pr = 0; pr < 2; pr ++) {
		if ((pr == 0) && (i == ii) && (j == jj)) {
			for (int jjj = j - 1; jjj <= j + 1; jjj++) {
				if (jjj == j)
					omega = 0.75;
				else
					omega = 0.125;
				lambda += omega * (multiplier(I, p, pr, ii + 1, jjj, kk, object_multiplier) - multiplier(I, p, pr, ii - 1, jjj, kk, object_multiplier));
			}
		} else if ((pr == 1) && (i == ii) && (j == jj)) {
			for (int iii = i - 1; iii <= i + 1; iii++) {
				if (iii == i)
					omega = 0.75;
				else
					omega = 0.125;
				lambda += omega * (multiplier(I, p, pr, iii, jj + 1, kk, object_multiplier) - multiplier(I, p, pr, iii, jj - 1, kk, object_multiplier));
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
			lambda += sign * omega * multiplier(I, p, pr, ii, jj, kk, object_multiplier);
		}
	}
	return lambda;
}

double delta(in *I, int p, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50])
{
	if ((i == ii) && (j == jj) && (k == kk))
		return 0;
	else
		return max3(0, - lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier), - lambda2D(I, p, ii, jj, kk, i, j, k, object_multiplier));
}

double alpha(in *I, int p, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_operand[50])
{
	if ((i == ii) && (j == jj) && (k == kk))
		return 0;
	double Q1 = 0, Q2 = 0, P1 = 0, P2 = 0;
	int kkk = k;
	for (int iii = i - 1; iii < i + 2; iii++) {
		for (int jjj = j - 1; jjj < j + 2; jjj++) {
			if (!((iii == 0) && (jjj == 0))) {
				Q1 += max2(0, lambda2D(I, p, i, j, k, iii, jjj, kkk, object_multiplier)) * max2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, i, j, k, object_operand), 0);
				Q2 += max2(0, lambda2D(I, p, i, j, k, iii, jjj, kkk, object_multiplier)) * min2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, i, j, k, object_operand), 0);
				P1 += min2(0, lambda2D(I, p, i, j, k, iii, jjj, kkk, object_multiplier)) * min2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, i, j, k, object_operand), 0);
				P2 += min2(0, lambda2D(I, p, i, j, k, iii, jjj, kkk, object_multiplier)) * max2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, i, j, k, object_operand), 0);
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
	if (lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier) <= lambda2D(I, p, ii, jj, kk, i, j, k, object_multiplier)) {
		if (operand(I, p, i, j, k, object_operand) >= operand(I, p, ii, jj, kk, object_operand))
			return min2(R1 * delta(I, p, i, j, k, ii, jj, kk, object_multiplier), lambda2D(I, p, ii, jj, kk, i, j, k, object_multiplier) + delta(I, p, ii, jj, kk, i, j, k, object_multiplier));
		else
			return min2(R2 * delta(I, p, i, j, k, ii, jj, kk, object_multiplier), lambda2D(I, p, ii, jj, kk, i, j, k, object_multiplier) + delta(I, p, ii, jj, kk, i, j, k, object_multiplier));
	} else {
		Q1 = Q2 = P1 = P2 = 0;
		int kkk = kk;
		for (int iii = ii - 1; iii < ii + 2; iii++) {
			for (int jjj = jj - 1; jjj < jj + 2; jjj++) {
				if (!((iii == 0) && (jjj == 0))) {
					Q1 += max2(0, lambda2D(I, p, ii, jj, kk, iii, jjj, kkk, object_multiplier)) * max2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, ii, jj, kk, object_operand), 0);
					Q2 += max2(0, lambda2D(I, p, ii, jj, kk, iii, jjj, kkk, object_multiplier)) * min2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, ii, jj, kk, object_operand), 0);
					P1 += min2(0, lambda2D(I, p, ii, jj, kk, iii, jjj, kkk, object_multiplier)) * min2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, ii, jj, kk, object_operand), 0);
					P2 += min2(0, lambda2D(I, p, ii, jj, kk, iii, jjj, kkk, object_multiplier)) * max2(operand(I, p, iii, jjj, kkk, object_operand) - operand(I, p, ii, jj, kk, object_operand), 0);
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
		if (operand(I, p, ii, jj, kk, object_operand) >= operand(I, p, i, j, k, object_operand))
			return min2(R1 * delta(I, p, ii, jj, kk, i, j, k, object_multiplier), lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier) + delta(I, p, i, j, k, ii, jj, kk, object_multiplier));
		else
			return min2(R2 * delta(I, p, ii, jj, kk, i, j, k, object_multiplier), lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier) + delta(I, p, i, j, k, ii, jj, kk, object_multiplier));
	}
}

double lambda_res(in *I, int p, int i, int j, int k, int ii, int jj, int kk, char object_multiplier[50], char object_operand[50], int antidiff)
{
	if ((i == ii) && (j == jj) && (k == kk))
		return lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier);
	else if (antidiff == 2)
		return lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier) + delta(I, p, i, j, k, ii, jj, kk, object_multiplier) - alpha(I, p, i, j, k, ii, jj, kk, object_multiplier, object_operand);
	else if (antidiff == 1)
		return lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier) + delta(I, p, i, j, k, ii, jj, kk, object_multiplier);
	else if (antidiff == 0)
		return lambda2D(I, p, i, j, k, ii, jj, kk, object_multiplier);
	return 0;
}

int div_func(in *I, int p, int i, int j, int k, char object_multiplier[50], char object_operand[50], int antidiff, double tetta)
{
	if (k != 0)
		return 0;
	double A_value, tmp;
	int kk = k;
	for (int ii = i - 1; ii < i + 2; ii++) {
		for (int jj = j - 1; jj < j + 2; jj++) {
			//if (!boundary_cell(I, ii, jj, kk)) {
			if (1) {
				tmp = lambda_res(I, p, i, j, k, ii, jj, kk, object_multiplier, object_operand, antidiff);
				A_value = tetta * tmp;
				WRITE_TO_A(p, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= (1 - tetta) * tmp * operand(I, p, i, j, k, object_operand);
				if (!((i == ii) && (j == jj) && (k == kk))) {
					A_value *= -1;
					WRITE_TO_A(p, ii, jj, kk, -1);
					I->B[A_IND(I, p, i, j, k)] += (1 - tetta) * tmp * operand(I, p, ii, jj, kk, object_operand);
				}
			}
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
	if (div_func(I, p, i, j, k, "density_avarage_velocity", "concentration", 1, 0.5)) return 1;
	return 0;
}

int SCAL_mass_inflow_rate_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	if ((p == 6) && (p == 7)) {
		A_value = - mass_inflow_rate_func(I, p, i, j, k) / saturation(I, p - 5, i, j, k);
		WRITE_TO_A(p, i, j, k, -1);
	} else if (p == 1) {
		A_value = - mass_inflow_rate_func(I, p, i, j, k) / concentration(I, p, i, j, k);
		WRITE_TO_A(p, i, j, k, -1);
	} else {
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, p, i, j, k);
	}
	return 0;
}

int LAPL_coef_pressure_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))) {
			A_value = (Darsi_M_coef(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) + Darsi_M_coef(I, i, j, k)) / (2 * I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			A_value = - A_value;
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		}
		if (!(boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			A_value = (Darsi_M_coef(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) + Darsi_M_coef(I, i, j, k)) / (2 * I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			A_value = - A_value;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

int DIV_density_average_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 5) && (p != 6) && (p != 7)) {
		printf("Incorrect index of DIV_density_average_velocity\n");
		return 1;
	}
	if (div_func(I, p, i, j, k, "density_avarage_velocity_divide_by_saturation", "saturation", 1, 0.5)) return 1;
	/*for (int pr = 0; pr < 3; pr++) {
		int ind_pr[3];
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		I->B[A_IND(I, p, i, j, k)] -= (density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
			density_t(I, p - 5, i, j, k) * avarage_velocity(I, p - 5, pr, i, j, k)) / (I->dx[pr] * I->porousness);
	}*/
	return 0;
}

int DIV_density_saturation_internal_energy_avarage_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, B_value;
	int pr, ind_pr[3], pp, ind_pp[3];
	double lambda1, lambda2, delta, alpha, Q1, Q2, P1, P2, R1, R2, tetta = 0.5;
	double multiplier(int ii, int jj, int kk)
	{
		return density_t(I, 2, ii, jj, kk) * avarage_velocity(I, 2, pr, ii, jj, kk) / (2 * I->dx[pr] * saturation(I, p - 5, ii, jj, kk));
	}
	double operand(int ii, int jj, int kk)
	{
		return saturation(I, p - 5, ii, jj, kk);
	}
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;

		A_value = B_value = 0;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			for (pp = 0; pp < 2; pp++) {
				A_value += density_t(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					I->specific_heat[pp] * I->adiabatic_exponent[pp] *
					avarage_velocity(I, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
				B_value += density_t(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					(- I->specific_heat[pp] * I->adiabatic_exponent[pp] * I->tempetarure_for_calculation_internal_energy + I->initial_enthalpy[pp]) *
					avarage_velocity(I, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			}
			for (pp = 0; pp < 4; pp++) {
				A_value += density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] *
					concentration(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
				B_value += density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					(- I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] * I->tempetarure_for_calculation_internal_energy + I->initial_enthalpy[pp + 2]) *
					concentration(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			}
			A_value /= (2 * I->dx[pr]);
			//printf("%d:\t%d\t%d\t", pr, i, k);
			//printf("%f\t", A_value);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			I->B[A_IND(I, p, i, j, k)] -= B_value;
			A_value = B_value = 0;
			for (pp = 0; pp < 2; pp++) {
				A_value -= density_t(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					I->specific_heat[pp] * I->adiabatic_exponent[pp] *
					avarage_velocity(I, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
				B_value -= density_t(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					(- I->specific_heat[pp] * I->adiabatic_exponent[pp] * I->tempetarure_for_calculation_internal_energy + I->initial_enthalpy[pp]) *
					avarage_velocity(I, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			}
			for (pp = 0; pp < 4; pp++) {
				A_value -= density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] *
					concentration(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
				B_value -= density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					(- I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] * I->tempetarure_for_calculation_internal_energy + I->initial_enthalpy[pp + 2]) *
					concentration(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			}
			A_value /= (2 * I->dx[pr]);
			//printf("%f\n", A_value);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			I->B[A_IND(I, p, i, j, k)] -= B_value;
		}
	}
	return 0;
}

int DIV_density_internal_energy_avarage_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if (p != 8) {
		printf("Incorrect index of DIV_density_internal_energy_average_velocity\n");
		return 1;
	}
	if (div_func(I, p, i, j, k, "density_internal_energy_avarage_velocity", "temperature_flow", 1, 0.5)) return 1;
	return 0;
}

int DIV_heat_influx_vector_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	for (int pr = 0; pr < 3; pr++) {
		int ind_pr[3];
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (! boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2])) {
			double A_value = I->porousness * (thermal_conductivity(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) + thermal_conductivity(I, i, j, k)) / (2 * I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			A_value *= -1;
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		}
		if (! boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) {
			double A_value = I->porousness * (thermal_conductivity(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) + thermal_conductivity(I, i, j, k)) / (2 * I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			A_value *= -1;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

int LAPL_thermal_conductivity_coef_saturation_temperature_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		A_value = I->porousness * (
			I->thermal_conductivity_coef[0] * saturation(I, 0, i, j, k) +
			I->thermal_conductivity_coef[1] * saturation(I, 1, i, j, k) +
			I->thermal_conductivity_coef[2] * saturation(I, 2, i, j, k)) / (I->dx[pr] * I->dx[pr]);
		WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		A_value *= -2;
		WRITE_TO_A(p, i, j, k, -1);
	}
	return 0;
}

int SCAL_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	I->B[A_IND(I, p, i, j, k)] += I->heat_transfer_coef * temperature_environment(I, i, j, k);
	A_value = I->heat_transfer_coef;
	WRITE_TO_A(p, i, j, k, -1);
	return 0;
}

int SCAL_chemical_reaction_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	//I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, 1, i, j, k) * enthalpy_flow(I, i, j, k);
	double A_value;
	int a;
//	A_value = - mass_inflow_rate_func(I, 1, i, j, k) *
//		(I->specific_heat[0] - I->specific_heat[1] -
//		 I->specific_heat[2] * concentration(I, 0, i, j, k) - I->specific_heat[3] * concentration(I, 1, i, j, k) +
//		 I->specific_heat[4] * concentration(I, 2, i, j, k) + I->specific_heat[5] * concentration(I, 3, i, j, k));
//	A_value = - mass_inflow_rate_func(I, 1, i, j, k) *
//		(4,5 * I->specific_heat[0] * saturation(I, 0, i, j, k) -
//		 I->specific_heat[1] * saturation(I, 1, i, j, k) -
//		 12.5 * I->specific_heat[3] * I->adiabatic_exponent[3] * saturation(I, 2, i, j, k) * concentration(I, 1, i, j, k) +
//		 8 * I->specific_heat[4] * I->adiabatic_exponent[4] * saturation(I, 2, i, j, k) * concentration(I, 2, i, j, k) +
//		 4.5 * I->specific_heat[5] * I->adiabatic_exponent[5] * saturation(I, 2, i, j, k) * concentration(I, 3, i, j, k));
//	WRITE_TO_A(p, i, j, k, -1);
//	I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func(I, 1, i, j, k) * (I->tempetarure_for_calculation_internal_energy *
//		(- I->specific_heat[0] + I->specific_heat[1] +
//		 I->specific_heat[2] * concentration(I, 0, i, j, k) + I->specific_heat[3] * concentration(I, 1, i, j, k) -
//		 I->specific_heat[4] * concentration(I, 2, i, j, k) - I->specific_heat[5] * concentration(I, 3, i, j, k)) +
//		I->initial_enthalpy[0] - I->initial_enthalpy[1] -
//		I->initial_enthalpy[2] * concentration(I, 0, i, j, k) - I->initial_enthalpy[3] * concentration(I, 1, i, j, k) +
//		I->initial_enthalpy[4] * concentration(I, 2, i, j, k) + I->initial_enthalpy[5] * concentration(I, 3, i, j, k));
//	I->B[A_IND(I, p, i, j, k)] += A_value * temperature_flow(I, i, j, k);
//	A_value = chemical_reaction_heat_flow_derivative_by_temperature(I, i, j, k);
//	WRITE_TO_A(p, i, j, k, -1);
//	I->B[A_IND(I, p, i, j, k)] += chemical_reaction_heat_flow_derivative_by_temperature(I, i, j, k) * temperature_flow(I, i, j, k) - chemical_reaction_heat_flow(I, i, j, k);
	/*	
 	A_value = 0;
	for (a = 0; a < 4; a++) {
		A_value += mass_inflow_rate_func_derivative_by_temperature(I, a, i, j, k) * I->initial_enthalpy[a + 2];
	}
	for (a = 0; a < 2; a++) {
		A_value += mass_inflow_rate_func_derivative_by_temperature(I, a + 5, i, j, k) * I->initial_enthalpy[a];
	}
	WRITE_TO_A(p, i, j, k, -1);
	*/
	for (a = 0; a < 4; a++) {
		I->B[A_IND(I, p, i, j, k)] -= mass_inflow_rate_func(I, a, i, j, k) * I->initial_enthalpy[a + 2];
	}
	for (a = 0; a < 2; a++) {
		I->B[A_IND(I, p, i, j, k)] -= mass_inflow_rate_func(I, a + 5, i, j, k) * I->initial_enthalpy[a];
	}
	/*
	for (a = 0; a < 4; a++) {
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func_derivative_by_temperature(I, a, i, j, k) * I->initial_enthalpy[a + 2] * temperature_flow(I, i, j, k);
	}
	for (a = 0; a < 2; a++) {
		I->B[A_IND(I, p, i, j, k)] += mass_inflow_rate_func_derivative_by_temperature(I, a + 5, i, j, k) * I->initial_enthalpy[a] * temperature_flow(I, i, j, k);
	}
	*/
	return 0;
}

int DIV_heat_influx_vector_environment_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	if (LAPL(p, i, j, k, thermal_conductivity_coef_saturation_temperature_environment, backward_euler, second, separated, FDM, termogas)) return 1;
	return 0;
}

int LAPL_thermal_conductivity_coef_saturation_temperature_environment_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		if (pr != 1) {
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			if (! boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2])) {
				A_value = (1 - I->porousness) * I->thermal_conductivity_coef[3] / (I->dx[pr] * I->dx[pr]);
				WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
			}
			if (! boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) {
				A_value = (1 - I->porousness) * I->thermal_conductivity_coef[3] / (I->dx[pr] * I->dx[pr]);
				WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
			}
		}
	}
	return 0;
}

int SCAL_minus_heat_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = I->heat_transfer_coef;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->heat_transfer_coef * temperature_flow(I, i, j, k);
	return 0;
}

