#include "init_data.h"
#include "utils.h"
#include "x_backward_euler_second_combined_FDM.h"
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

int DIV_grad_pressure_backward_euler_second_combined_FDM(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		if (pr != 1) {
			ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
			ind_pr[pr] = 1;
			A_value = 1 / (I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
			A_value = - 2 / (I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i, j, k, -1);
			A_value = 1 / (I->dx[pr] * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

