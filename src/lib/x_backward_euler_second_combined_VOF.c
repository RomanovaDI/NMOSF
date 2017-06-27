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

int DIV_grad_pressure_backward_euler_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int s;
	double A_value;
	for (s = 0; s < 6; s++) {
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 2.) * I->normal[NORMAL_IND(I, 0, i, j, k, s)] / (2. * I->dx[0]);
		WRITE_TO_A(4, i + 1, j, k, -1);
		WRITE_TO_A(4, i + 1, j, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i - 1, j, k, -1);
		WRITE_TO_A(4, i - 1, j, k, s);
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 2.) * I->normal[NORMAL_IND(I, 1, i, j, k, s)] / (2. * I->dx[1]);
		WRITE_TO_A(4, i, j + 1, k, -1);
		WRITE_TO_A(4, i, j + 1, k, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j - 1, k, -1);
		WRITE_TO_A(4, i, j - 1, k, s);
		A_value = (I->area[AREA_IND(I, i, j, k, s)] / 2.) * I->normal[NORMAL_IND(I, 2, i, j, k, s)] / (2. * I->dx[2]);
		WRITE_TO_A(4, i, j, k + 1, -1);
		WRITE_TO_A(4, i, j, k + 1, s);
		A_value = - A_value;
		WRITE_TO_A(4, i, j, k - 1, -1);
		WRITE_TO_A(4, i, j, k - 1, s);
	}
	return 0;
}

