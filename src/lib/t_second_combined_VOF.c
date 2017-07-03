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
	//I->B[A_IND(I, p, i, j, k)] += (
	//	density(I, i, j, k) * velocity(I, p, i, j, k) +
	//	density(I, i - 1, j, k) * velocity(I, p, i - 1, j, k) +
	//	density(I, i + 1, j, k) * velocity(I, p, i + 1, j, k) +
	//	density(I, i, j - 1, k) * velocity(I, p, i, j - 1, k) +
	//	density(I, i, j + 1, k) * velocity(I, p, i, j + 1, k) +
	//	density(I, i, j, k - 1) * velocity(I, p, i, j, k - 1) +
	//	density(I, i, j, k + 1) * velocity(I, p, i, j, k + 1)) / (7. * I->dt);
	I->B[A_IND(I, p, i, j, k)] += density(I, i, j, k) * (
		//velocity(I, p, i, j, k) +
		velocity(I, p, i - 1, j, k) +
		velocity(I, p, i + 1, j, k) +
		velocity(I, p, i, j - 1, k) +
		velocity(I, p, i, j + 1, k) +
		velocity(I, p, i, j, k - 1) +
		velocity(I, p, i, j, k + 1)) / (6. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	density(I, i, j, k) * velocity(I, p, i, j, k) +
	//	density(I, i - 1, j, k) * velocity(I, p, i - 1, j, k) +
	//	density(I, i + 1, j, k) * velocity(I, p, i + 1, j, k) +
	//	density(I, i, j, k - 1) * velocity(I, p, i, j, k - 1) +
	//	density(I, i, j, k + 1) * velocity(I, p, i, j, k + 1)) / (4. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	density(I, i, j, k) * velocity(I, p, i, j, k) +
	//	density(I, i - 1, j, k) * velocity(I, p, i - 1, j, k) +
	//	density(I, i + 1, j, k) * velocity(I, p, i + 1, j, k) +
	//	density(I, i, j - 1, k) * velocity(I, p, i, j - 1, k) +
	//	density(I, i, j + 1, k) * velocity(I, p, i, j + 1, k) +
	//	density(I, i, j, k - 1) * velocity(I, p, i, j, k - 1) +
	//	density(I, i, j, k + 1) * velocity(I, p, i, j, k + 1)) / (6. * I->dt);
	return 0;
}

int DDT_snow_volume_fraction_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = 1 / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	I->B[A_IND(I, p, i, j, k)] += (
		//phase_fraction(I, i, j, k) +
		phase_fraction(I, i - 1, j, k) +
        phase_fraction(I, i + 1, j, k) +
        phase_fraction(I, i, j - 1, k) +
        phase_fraction(I, i, j + 1, k) +
        phase_fraction(I, i, j, k - 1) +
        phase_fraction(I, i, j, k + 1)) / (6. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	phase_fraction(I, i, j, k) +
	//	phase_fraction(I, i - 1, j, k) +
    //    phase_fraction(I, i + 1, j, k) +
    //    phase_fraction(I, i, j, k - 1) +
    //    phase_fraction(I, i, j, k + 1)) / (4. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	phase_fraction(I, i, j, k) +
	//	phase_fraction(I, i - 1, j, k) +
    //    phase_fraction(I, i + 1, j, k) +
    //    phase_fraction(I, i, j - 1, k) +
    //    phase_fraction(I, i, j + 1, k) +
    //    phase_fraction(I, i, j, k - 1) +
    //    phase_fraction(I, i, j, k + 1)) / (6. * I->dt);
	return 0;
}

int DDT_pressure_cont_second_combined_VOF(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, prerho, c1 = 1, c2 = 0;
	prerho = (I->density_snow / I->pressure_atmosphere) * (c1 + c2 * (pressure(I, i, j, k) - I->pressure_atmosphere));
	//prerho = I->density_snow / I->pressure_atmosphere;
	A_value = prerho / I->dt;
	WRITE_TO_A(p, i, j, k, -1);
	//I->B[A_IND(I, p, i, j, k)] += prerho * (
	//	pressure(I, i, j, k) +
	//	pressure(I, i - 1, j, k) +
    //    pressure(I, i + 1, j, k) +
    //    pressure(I, i, j - 1, k) +
    //    pressure(I, i, j + 1, k) +
    //    pressure(I, i, j, k - 1) +
    //    pressure(I, i, j, k + 1)) / (7. * I->dt);
	I->B[A_IND(I, p, i, j, k)] += prerho * (
	//	pressure(I, i, j, k) +
		pressure(I, i - 1, j, k) +
        pressure(I, i + 1, j, k) +
    //    pressure(I, i, j - 1, k) +
    //    pressure(I, i, j + 1, k) +
        pressure(I, i, j, k - 1) +
        pressure(I, i, j, k + 1)) / (4. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	phase_fraction(I, i, j, k) +
	//	phase_fraction(I, i - 1, j, k) +
    //    phase_fraction(I, i + 1, j, k) +
    //    phase_fraction(I, i, j, k - 1) +
    //    phase_fraction(I, i, j, k + 1)) / (4. * I->dt);
	//I->B[A_IND(I, p, i, j, k)] += (
	////	phase_fraction(I, i, j, k) +
	//	phase_fraction(I, i - 1, j, k) +
    //    phase_fraction(I, i + 1, j, k) +
    //    phase_fraction(I, i, j - 1, k) +
    //    phase_fraction(I, i, j + 1, k) +
    //    phase_fraction(I, i, j, k - 1) +
    //    phase_fraction(I, i, j, k + 1)) / (6. * I->dt);
	return 0;
}
