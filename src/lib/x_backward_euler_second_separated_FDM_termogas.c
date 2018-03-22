#include "init_data.h"
#include "utils.h"
#include "x_backward_euler_second_separated_FDM_termogas.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

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

double lambda(in *I, int p, int pr, int i, int j, int k, int ii, int jj, int kk, int iii, int jjj, int kkk)
{
	if ((i == ii) && (j == jj) && (k == kk))
		int omega = 0.75;
	else
		int omega = 0.125;
	if ()
}

int ddx_calc(in *I, int p, int pr, int i, int j, int k, double omega, double tetta, int antidiff, int ind) // i_n, j_n, k_n - coordinates of neighbor cell
{
	double A_value;
	int pp, ind_pp[3], ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	double lambda1, lambda2, delta, alpha, Q1, Q2, P1, P2, R1, R2, RR1, RR2;
	double multiplier(int ii, int jj, int kk)
	{
		switch (ind) {
			case 0:
				return density_t(I, 2, ii, jj, kk) * avarage_velocity(I, 2, pr, ii, jj, kk) / (2 * I->dx[pr]);
			case 1:
				return density_t(I, p - 5, ii, jj, kk) * avarage_velocity(I, p - 5, pr, ii, jj, kk) / (2 * I->dx[pr] * saturation(I, p - 5, ii, jj, kk));
		}
	}
	double operand(int ii, int jj, int kk)
	{
		switch (ind){
			case 0:
				return concentration(I, p, ii, jj, kk);
			case 1:
				return saturation(I, p - 5, ii, jj, kk);
		}
	}
	lambda1 = - omega * multiplier(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
	lambda2 = omega * multiplier(i, j, k);
	delta = max3(0, -lambda1, -lambda2);
	if (antidiff) {
		Q1 = Q2 = P1 = P2 = 0;
		for (pp = 0; pp < 2; pp++) {
			ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
			ind_pp[pp] = 1;
			Q1 += max2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * max2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
			Q1 += max2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * max2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
			Q2 += max2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * min2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
			Q2 += max2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * min2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
			P1 += min2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * min2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
			P1 += min2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * min2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
			P2 += min2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * max2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
			P2 += min2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * max2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
		}
		if (P1 == 0)
			Q1 = 0;
		else
			Q1 /= P1;
		if (P2 == 0)
			Q2 = 0;
		else
			Q2 /= P2;
		R1 = max2(0, min2(1, Q1));
		R2 = max2(0, min2(1, Q2));
		if (lambda1 <= lambda2) {
			if (operand(i, j, k) >= operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
				alpha = min2(R1 * delta, lambda2 + delta);
			else
				alpha = min2(R2 * delta, lambda2 + delta);
		} else {
			Q1 = Q2 = P1 = P2 = 0;
			for (pp = 0; pp < 2; pp++) {
				ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
				ind_pp[pp] = 1;
				Q1 += max2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * max2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				Q1 += max2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * max2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				Q2 += max2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * min2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				Q2 += max2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * min2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				P1 += min2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * min2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				P1 += min2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * min2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				P2 += min2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * max2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				P2 += min2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * max2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
			}
			if (P1 == 0)
				Q1 = 0;
			else
				Q1 /= P1;
			if (P2 == 0)
				Q2 = 0;
			else
				Q2 /= P2;
			RR1 = max2(0, min2(1, Q1));
			RR2 = max2(0, min2(1, Q2));
			if (operand(i, j, k) <= operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
				alpha = min2(RR1 * delta, lambda1 + delta);
			else
				alpha = min2(RR2 * delta, lambda1 + delta);
		}
	} else
		alpha = 0;
	A_value = -lambda1 - delta + alpha;
	A_value *= tetta;
	WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
	I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
	A_value *= -1;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);

	lambda1 = multiplier(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
	lambda2 = -multiplier(i, j, k);
	delta = max3(0, -lambda1, -lambda2);
	if (antidiff) {
		if (lambda1 <= lambda2) {
			if (operand(i, j, k) >= operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
				alpha = min2(R1 * delta, lambda2 + delta);
			else
				alpha = min2(R2 * delta, lambda2 + delta);
		} else {
			Q1 = Q2 = P1 = P2 = 0;
			for (pp = 0; pp < 2; pp++) {
				ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
				ind_pp[pp] = 1;
				Q1 += max2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * max2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				Q1 += max2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * max2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				Q2 += max2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * min2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				Q2 += max2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * min2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				P1 += min2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * min2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				P1 += min2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * min2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				P2 += min2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * max2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				P2 += min2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * max2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
			}
			if (P1 == 0)
				Q1 = 0;
			else
				Q1 /= P1;
			if (P2 == 0)
				Q2 = 0;
			else
				Q2 /= P2;
			RR1 = max2(0, min2(1, Q1));
			RR2 = max2(0, min2(1, Q2));
			if (operand(i, j, k) <= operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
				alpha = min2(RR1 * delta, lambda1 + delta);
			else
				alpha = min2(RR2 * delta, lambda1 + delta);
		}
	} else
		alpha = 0;
	A_value = -lambda1 - delta + alpha;
	A_value *= tetta;
	WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
	I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
	A_value *= -1;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);

	A_value = multiplier(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - multiplier(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
	A_value *= tetta;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);
	return 0;
}

int DIV_concentration_density_average_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double tmp, a;
	double A_value;
	int pr, ind_pr[3], pp, ind_pp[3];
	for (pr = 0; pr < 2; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (1) {
		//if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
/*
			A_value = density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
#if DEBUG
			printf("density_t(2, %d, %d, %d) = %f\t avarage_velocity(I, 2, %d, %d, %d, %d) = %f\n",\
				i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]),\
				pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]));
#endif
			A_value = - density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
#if DEBUG
			printf("density_t(2, %d, %d, %d) = %f\t avarage_velocity(I, 2, %d, %d, %d, %d) = %f\n",\
				i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]),\
				pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]));
#endif
*/
/*
			tmp = density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (2 * I->dx[pr]);
			tmp *= - density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
			if (tmp > 0)
				a = 2.0;
			else
				a = 1.0;
			A_value = 0.5 * density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (a * I->dx[pr]);
			if (A_value < 0) {
				WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * concentration(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * concentration(I, p, i, j, k);
			}
			A_value = - 0.5 * density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (a * I->dx[pr]);
			if (A_value < 0) {
				WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * concentration(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * concentration(I, p, i, j, k);
			}
			A_value = 0.5 * (density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
				density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * concentration(I, p, i, j, k);
*/
/*
			A_value = 1.0 / (I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			//WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value *= -1;
			//WRITE_TO_A(p, i, j, k, -1);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
*/

			if (ddx_calc(I, p, pr, i, j, k, 1, 0.5, 0, 0)) return 1;

/*
			if (pr == 0)
				pp = 1;
			else if (pr == 1)
				pp = 0;
			ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
			ind_pp[pp] = 1;
			tmp = 0.75;
			if (pr != 2) {
				if	(internal_cell(I, i + ind_pr[0] + ind_pp[0], j + ind_pr[1] + ind_pp[1], k + ind_pr[2] + ind_pp[2]) &&
					internal_cell(I, i - ind_pr[0] + ind_pp[0], j - ind_pr[1] + ind_pp[1], k - ind_pr[2] + ind_pp[2])) {
						if (ddx_calc(I, p, pr, i + ind_pp[0], j + ind_pp[1], k + ind_pp[2], 0.125, 0.5, 0, 0)) return 1;
				} else
					tmp += 0.125;
				if	(internal_cell(I, i + ind_pr[0] - ind_pp[0], j + ind_pr[1] - ind_pp[1], k + ind_pr[2] - ind_pp[2]) &&
					internal_cell(I, i - ind_pr[0] - ind_pp[0], j - ind_pr[1] - ind_pp[1], k - ind_pr[2] - ind_pp[2])) {
						if (ddx_calc(I, p, pr, i - ind_pp[0], j - ind_pp[1], k - ind_pp[2], 0.125, 0.5, 0, 0)) return 1;
				} else
					tmp += 0.125;
			}
			if (ddx_calc(I, p, pr, i, j, k, tmp, 0.5, 0, 0)) return 1;
*/
		}
	}
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
//		A_value = - Darsi_M_coef(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (I->dx[pr] * I->dx[pr]);
//		WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
//		A_value = - Darsi_M_coef(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (I->dx[pr] * I->dx[pr]);
//		WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
//		A_value *= -2;
//		WRITE_TO_A(p, i, j, k, -1);
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
	double A_value, tmp, a;
	int pr, ind_pr[3], pp, ind_pp[3];
	for (pr = 0; pr < 2; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (1) {
		//if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
/*
			I->B[A_IND(I, p, i, j, k)] -= (density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
				density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
*/
/*
			tmp = density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
				(2 * I->dx[pr] * saturation(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]));
			tmp *= - density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
				(2 * I->dx[pr] * saturation(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]));
			//tmp = density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
			//	(2 * I->dx[pr]);
			//tmp *= - density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
			//	(2 * I->dx[pr]);
			if (tmp > 0)
				a = 2.0;
			else
				a = 1.0;
			A_value = 0.5 * density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
				(a * I->dx[pr] * saturation(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]));
			//A_value = 0.5 * density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
			//	(a * I->dx[pr]);
			if (A_value < 0) {
				WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * saturation(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
				//I->B[A_IND(I, p, i, j, k)] -= A_value;
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * saturation(I, p - 5, i, j, k);
				//I->B[A_IND(I, p, i, j, k)] -= A_value;
			}
			A_value = - 0.5 * density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
				(a * I->dx[pr] * saturation(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]));
			//A_value = - 0.5 * density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
			//	(a * I->dx[pr]);
			if (A_value < 0) {
				WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * saturation(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
				//I->B[A_IND(I, p, i, j, k)] -= A_value;
				A_value *= -1;
				WRITE_TO_A(p, i, j, k, -1);
				I->B[A_IND(I, p, i, j, k)] -= A_value * saturation(I, p - 5, i, j, k);
				//I->B[A_IND(I, p, i, j, k)] -= A_value;
			}
			A_value = 0.5 * (density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
				saturation(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
				density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
				saturation(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) /
				(2 * I->dx[pr]);
			//A_value = 0.5 * (density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
			//	density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
			//	avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) /
			//	(2 * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * saturation(I, p - 5, i, j, k);
			//I->B[A_IND(I, p, i, j, k)] -= A_value;
*/
/*
			A_value = density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) /
				(saturation(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * 2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value = -	density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
				avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) /
				(saturation(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * 2 * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
*/
/*
			A_value = 1.0 / (2.0 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value *= -1;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
*/

			if (ddx_calc(I, p, pr, i, j, k, 1, 0.5, 0, 1)) return 1;

/*
			if (pr == 0)
				pp = 1;
			else if (pr == 1)
				pp = 0;
			ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
			ind_pp[pp] = 1;
			tmp = 0.75;
			if (pr != 2) {
				if	(internal_cell(I, i + ind_pr[0] + ind_pp[0], j + ind_pr[1] + ind_pp[1], k + ind_pr[2] + ind_pp[2]) &&
					internal_cell(I, i - ind_pr[0] + ind_pp[0], j - ind_pr[1] + ind_pp[1], k - ind_pr[2] + ind_pp[2])) {
						if (ddx_calc(I, p, pr, i + ind_pp[0], j + ind_pp[1], k + ind_pp[2], 0.125, 0.5, 0, 1)) return 1;
				} else
					tmp += 0.125;
				if	(internal_cell(I, i + ind_pr[0] - ind_pp[0], j + ind_pr[1] - ind_pp[1], k + ind_pr[2] - ind_pp[2]) &&
					internal_cell(I, i - ind_pr[0] - ind_pp[0], j - ind_pr[1] - ind_pp[1], k - ind_pr[2] - ind_pp[2])) {
						if (ddx_calc(I, p, pr, i - ind_pp[0], j - ind_pp[1], k - ind_pp[2], 0.125, 0.5, 0, 1)) return 1;
				} else
					tmp += 0.125;
			}
			if (ddx_calc(I, p, pr, i, j, k, tmp, 0.5, 0, 1)) return 1;
*/
		}
	}
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

/*
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			lambda1 = - multiplier(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			lambda2 = multiplier(i, j, k);
			delta = max3(0, -lambda1, -lambda2);
			Q1 = Q2 = P1 = P2 = 0;
			for (pp = 0; pp < 3; pp++) {
				ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
				ind_pp[pp] = 1;
				Q1 += max2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * max2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
				Q1 += max2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * max2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
				Q2 += max2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * min2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
				Q2 += max2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * min2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
				P1 += min2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * min2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
				P1 += min2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * min2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
				P2 += min2(0, -multiplier(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2])) * max2(operand(i + ind_pp[0], j + ind_pp[1], k + ind_pp[2]) - operand(i, j, k), 0);
				P2 += min2(0,  multiplier(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2])) * max2(operand(i - ind_pp[0], j - ind_pp[1], k - ind_pp[2]) - operand(i, j, k), 0);
			}
			if (P1 == 0)
				Q1 = 0;
			else
				Q1 /= P1;
			if (P2 == 0)
				Q2 = 0;
			else
				Q2 /= P2;
			R1 = max2(0, min2(1, Q1));
			R2 = max2(0, min2(1, Q2));
			if (lambda1 <= lambda2) {
				if (operand(i, j, k) >= operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
					alpha = min2(R1 * delta, lambda2 + delta);
				else
					alpha = min2(R2 * delta, lambda2 + delta);
			} else {
				Q1 = Q2 = P1 = P2 = 0;
				for (pp = 0; pp < 3; pp++) {
					ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
					ind_pp[pp] = 1;
					Q1 += max2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * max2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					Q1 += max2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * max2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					Q2 += max2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * min2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					Q2 += max2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * min2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					P1 += min2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * min2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					P1 += min2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * min2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					P2 += min2(0, -multiplier(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2])) * max2(operand(i + ind_pp[0] + ind_pr[0], j + ind_pp[1] + ind_pr[1], k + ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
					P2 += min2(0,  multiplier(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2])) * max2(operand(i - ind_pp[0] + ind_pr[0], j - ind_pp[1] + ind_pr[1], k - ind_pp[2] + ind_pr[2]) - operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]), 0);
				}
				if (P1 == 0)
					Q1 = 0;
				else
					Q1 /= P1;
				if (P2 == 0)
					Q2 = 0;
				else
					Q2 /= P2;
				R1 = max2(0, min2(1, Q1));
				R2 = max2(0, min2(1, Q2));
				if (operand(i, j, k) <= operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
					alpha = min2(R1 * delta, lambda1 + delta);
				else
					alpha = min2(R2 * delta, lambda1 + delta);
			}
			A_value = -lambda1 - delta + alpha;
			A_value *= tetta;
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			A_value *= -1;
			WRITE_TO_A(p, i, j, k, -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);
			lambda1 = multiplier(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			lambda2 = -multiplier(i, j, k);
			delta = max3(0, -lambda1, -lambda2);
			if (lambda1 <= lambda2) {
				if (operand(i, j, k) >= operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
					alpha = min2(R1 * delta, lambda2 + delta);
				else
					alpha = min2(R2 * delta, lambda2 + delta);
			} else {
				Q1 = Q2 = P1 = P2 = 0;
				for (pp = 0; pp < 3; pp++) {
					ind_pp[0] = ind_pp[1] = ind_pp[2] = 0;
					ind_pp[pp] = 1;
					Q1 += max2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * max2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					Q1 += max2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * max2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					Q2 += max2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * min2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					Q2 += max2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * min2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					P1 += min2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * min2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					P1 += min2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * min2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					P2 += min2(0, -multiplier(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2])) * max2(operand(i + ind_pp[0] - ind_pr[0], j + ind_pp[1] - ind_pr[1], k + ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
					P2 += min2(0,  multiplier(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2])) * max2(operand(i - ind_pp[0] - ind_pr[0], j - ind_pp[1] - ind_pr[1], k - ind_pp[2] - ind_pr[2]) - operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]), 0);
				}
				if (P1 == 0)
					Q1 = 0;
				else
					Q1 /= P1;
				if (P2 == 0)
					Q2 = 0;
				else
					Q2 /= P2;
				R1 = max2(0, min2(1, Q1));
				R2 = max2(0, min2(1, Q2));
				if (operand(i, j, k) <= operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
					alpha = min2(R1 * delta, lambda1 + delta);
				else
					alpha = min2(R2 * delta, lambda1 + delta);
			}
			A_value = -lambda1 - delta + alpha;
			A_value *= tetta;
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			A_value *= -1;
			WRITE_TO_A(p, i, j, k, -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);
			A_value = multiplier(i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - multiplier(i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			A_value *= tetta;
			WRITE_TO_A(p, i, j, k, -1);
			I->B[A_IND(I, p, i, j, k)] -= A_value * operand(i, j, k);
*/
		}
	}
	return 0;
}

int DIV_density_internal_energy_avarage_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pp, pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		A_value = 0;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			for (pp = 0; pp < 2; pp++) {
				A_value += density_t(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					I->specific_heat[pp] * I->adiabatic_exponent[pp] *
					avarage_velocity(I, pp, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			}
			for (pp = 0; pp < 4; pp++) {
				A_value += density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] *
					concentration(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) *
					avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]);
			}
			A_value /= (2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value = 0;
			for (pp = 0; pp < 2; pp++) {
				A_value -= density_t(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					I->specific_heat[pp] * I->adiabatic_exponent[pp] *
					avarage_velocity(I, pp, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			}
			for (pp = 0; pp < 4; pp++) {
				A_value -= density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] *
					concentration(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) *
					avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]);
			}
			A_value /= (2 * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

int DIV_heat_influx_vector_flow_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, x;
	int pr, ind_pr[3], pp;
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (! boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2])) {
			A_value = - I->porousness;
			x = 0;
			for (pp = 0; pp < 3; pp++)
				x += I->thermal_conductivity_coef[pp] * (saturation(I, pp, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) + saturation(I, pp, i, j, k)) / 2;
			A_value *= x / (I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value *= -1;
			WRITE_TO_A(p, i, j, k, -1);
		}
		if (! boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) {
			A_value = - I->porousness;
			x = 0;
			for (pp = 0; pp < 3; pp++)
				x += I->thermal_conductivity_coef[pp] * (saturation(I, pp, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) + saturation(I, pp, i, j, k)) / 2;
			A_value *= x / (I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
			A_value *= -1;
			WRITE_TO_A(p, i, j, k, -1);
		}
	}
//	if (LAPL(p, i, j, k, thermal_conductivity_coef_saturation_temperature_flow, backward_euler, second, separated, FDM, termogas)) return 1;
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

