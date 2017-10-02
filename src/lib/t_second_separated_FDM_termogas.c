#include "init_data.h"
#include "utils.h"
#include "t_second_separated_FDM_termogas.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WRITE_TO_A(pr_m, i_m, j_m, k_m, s_m)\
	if (write_to_A_csr(I, p, i, j, k, pr_m, i_m, j_m, k_m, s_m, A_value)) { \
		printf("Error writing matrix A\n"); \
		return 1; \
	}

int DDT_concentration_density_saturation_porousness_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = density_t(I, p, i, j, k) * saturation(I, 2, i, j, k) * I->porousness / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
//	I->B[A_IND(I, p, i, j, k)] += A_value * (
//		concentration(I, p, i + 1, j, k) +
//		concentration(I, p, i - 1, j, k) +
//		concentration(I, p, i, j + 1, k) +
//		concentration(I, p, i, j - 1, k) +
//		concentration(I, p, i, j, k + 1) +
//		concentration(I, p, i, j, k - 1)) / 6;
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		concentration(I, p, i + 1, j, k) +
		concentration(I, p, i - 1, j, k) +
		concentration(I, p, i, j + 1, k) +
		concentration(I, p, i, j - 1, k)) / 4;
	return 0;
}

int DDT_coef_pressure_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = Darsi_A_coef(I, i, j, k) / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		pressure(I, i + 1, j, k) + pressure(I, i - 1, j, k) +
		pressure(I, i, j + 1, k) + pressure(I, i, j - 1, k) +
		pressure(I, i, j, k + 1) + pressure(I, i, j, k - 1)) / 6;
	return 0;
}

int DDT_density_saturation_porousness_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = density_t(I, p - 4, i, j, k) * I->porousness / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		saturation(I, p - 4, i + 1, j, k) + saturation(I, p - 4, i - 1, j, k) +
		saturation(I, p - 4, i, j + 1, k) + saturation(I, p - 4, i, j - 1, k) +
		saturation(I, p - 4, i, j, k + 1) + saturation(I, p - 4, i, j, k - 1)) / 6;
	return 0;
}

int DDT_porousness_density_energy_flow_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value = 0;
	int pp;
	for (pp = 0; pp < 2; pp++)
		A_value += I->porousness * density_t(I, pp, i, j, k) * saturation(I, pp, i, j, k) * I->specific_heat[pp];
	for (pp = 0; pp < 4; pp++)
		A_value += I->porousness * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[pp + 2] * concentration(I, pp, i, j, k);
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		temperature_flow(I, i + 1, j, k) + temperature_flow(I, i - 1, j, k) +
		temperature_flow(I, i, j + 1, k) + temperature_flow(I, i, j - 1, k) +
		temperature_flow(I, i, j, k + 1) + temperature_flow(I, i, j, k - 1)) / 6;
	return 0;
}

int DDT_porousness_density_energy_environment_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = (1 - I->porousness) * I->density_environment * I->specific_heat[6];
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += A_value * (
		temperature_environment(I, i + 1, j, k) + temperature_environment(I, i - 1, j, k) +
		temperature_environment(I, i, j + 1, k) + temperature_environment(I, i, j - 1, k) +
		temperature_environment(I, i, j, k + 1) + temperature_environment(I, i, j, k - 1)) / 6;
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

