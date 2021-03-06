#include "init_data.h"
#include "utils.h"
#include "t_second_separated_FDM_termogas.h"
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DDT_concentration_density_saturation_porousness_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 0) && (p != 1) && (p != 2) && (p != 3)) {
		printf("Incorrect index of DDT_concentration_density_saturation_porousness\n");
		return 1;
	}
	double A_value;
	//A_value = density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->porousness / I->dt;
	A_value = 1.0 / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * concentration(I, p, i, j, k);
	return 0;
}

int DDT_arithmetic_mean_of_neighboring_cells_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value = 1;
	WRITE_TO_A(p, i, j, k, -1);
	A_value = - (6. / 13.);
	if (internal_cell(I, i + 1, j, k))
		WRITE_TO_A(p, i + 1, j, k, -1);
	if (internal_cell(I, i - 1, j, k))
		WRITE_TO_A(p, i - 1, j, k, -1);
	if (internal_cell(I, i, j + 1, k))
		WRITE_TO_A(p, i, j + 1, k, -1);
	if (internal_cell(I, i, j - 1, k))
		WRITE_TO_A(p, i, j - 1, k, -1);
	A_value = - (1. / 13.);
	if (internal_cell(I, i + 1, j + 1, k))
		WRITE_TO_A(p, i + 1, j + 1, k, -1);
	if (internal_cell(I, i + 1, j - 1, k))
		WRITE_TO_A(p, i + 1, j - 1, k, -1);
	if (internal_cell(I, i - 1, j + 1, k))
		WRITE_TO_A(p, i - 1, j + 1, k, -1);
	if (internal_cell(I, i - 1, j - 1, k))
		WRITE_TO_A(p, i - 1, j - 1, k, -1);
	return 0;
}

int DDT_coef_pressure_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = Darsi_A_coef(I, i, j, k) / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * pressure(I, i, j, k);
	return 0;
}

int DDT_density_saturation_porousness_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 5) && (p != 6) && (p != 7)) {
		printf("Incorrect index of DDT_density_saturation_porousness\n");
		return 1;
	}
	double A_value;
	A_value = density_t(I, p - 5, i, j, k) * I->porousness[POR_IND(I, i, j, k)] / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (saturation(I, p - 5, i, j, k) - I->residual_saturation[p - 5]);
	return 0;
}

int DDT_saturation_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	if ((p != 5) && (p != 6) && (p != 7)) {
		printf("Incorrect index of DDT_density_saturation_porousness\n");
		return 1;
	}
	double A_value;
	A_value = 1.0 / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (saturation(I, p - 5, i, j, k) - I->residual_saturation[p - 5]);
	return 0;
}

int DDT_porousness_density_energy_flow_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = 1.0 / I->dt;
	/*int pp;
	for (pp = 0; pp < 1; pp++)
		A_value += I->porousness * density_t(I, pp, i, j, k) * saturation(I, pp, i, j, k) * I->specific_heat[pp] / I->dt;
	for (pp = 0; pp < 0; pp++)
		A_value += I->porousness * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[pp + 2] * concentration(I, pp, i, j, k) / I->dt;*/
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (temperature_flow(I, i, j, k) - I->initial_temperature);
	return 0;
}

int DDT_porousness_density_energy_environment_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	//A_value = (1 - I->porousness) * I->density_environment * I->specific_heat[6] / I->dt;
	A_value = 1.0 / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (temperature_environment(I, i, j, k) - I->initial_temperature);
	return 0;
}

int DDT_identity_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value = 1;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += I->B_prev[B_IND(I, p, i, j, k)];
	return 0;
}

