#include "init_data.h"
#include "utils.h"
#include "t_second_combined_VOF.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DDT_density_velocity_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = density(I, i, j, k) / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += density(I, i, j, k) * velocity(I, p, i, j, k) / I->dt;
	return 0;
}

int DDT_snow_volume_fraction_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = 1 / I->dt;
	WRITE_TO_A(3, i, j, k, -1);
	I->B[A_IND(I, 3, i, j, k)] += phase_fraction(I, i, j, k) / I->dt;
	return 0;
}
