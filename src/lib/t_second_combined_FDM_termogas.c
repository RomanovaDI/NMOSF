#include "init_data.h"
#include "utils.h"
#include "t_second_combined_FDM_termogas.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DDT_concentration_density_saturation_porousness_second_combined_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = density_t(I, p, i, j, k) * saturation(I, 2, i, j, k) * I->porousness / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		concentration(I, p, i + 1, j, k) +
		concentration(I, p, i - 1, j, k) +
		concentration(I, p, i, j + 1, k) +
		concentration(I, p, i, j - 1, k) +
		concentration(I, p, i, j, k + 1) +
		concentration(I, p, i, j, k - 1)) / (6 * I->dt);
	return 0;
}
