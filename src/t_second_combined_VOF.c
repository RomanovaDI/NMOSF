#include "utils.h"
#include "t_second_combined_VOF.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DDT_density_velocity_second_combined_VOF(int p, int i, int j, int k);
int DDT_density_snow_volume_fraction_second_combined_VOF(int p, int i, int j, int k);

int DDT_density_velocity_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	A_value = density(i, j, k) / dt;
	WRITE_TO_A(p, i, j, k, -1);
	B[A_IND(p, i, j, k)] += density(i, j, k) * velocity(p, i, j, k) / dt;
	return 0;
}

int DDT_density_snow_volume_fraction_second_combined_VOF(int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(i, j, k)) return 1;
	double A_value;
	A_value = 1 / dt;
	WRITE_TO_A(3, i, j, k, -1);
	B[A_IND(3, i, j, k)] += phase_fraction(i, j, k) / dt;
	return 0;
}

