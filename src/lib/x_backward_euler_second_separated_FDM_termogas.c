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

int DIV_concentration_density_average_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	int pr, ind_pr[3];
	double A_value;
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			A_value = density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (2 * I->dx[pr]);
			WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
//			printf("density_t(2, %d, %d, %d) = %f\t avarage_velocity(I, 2, %d, %d, %d, %d) = %f\n",\
//				i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], density_t(I, 2, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]),\
//				pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], avarage_velocity(I, 2, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]));
			A_value = - density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
			WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
//			printf("density_t(2, %d, %d, %d) = %f\t avarage_velocity(I, 2, %d, %d, %d, %d) = %f\n",\
//				i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], density_t(I, 2, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]),\
//				pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], avarage_velocity(I, 2, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]));
		}
	}
	return 0;
}

int SCAL_mass_inflow_rate_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value;
	A_value = - mass_inflow_rate_func(I, p, i, j, k);
	WRITE_TO_A(p, i, j, k, -1);
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
	double A_value;
	int pr, ind_pr[3];
	for (pr = 0; pr < 3; pr++) {
		ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
		ind_pr[pr] = 1;
		if (!(boundary_cell(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) || boundary_cell(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))) {
			I->B[A_IND(I, p, i, j, k)] -= (density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) -
				density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
		//	A_value = density_t(I, p - 5, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) * avarage_velocity(I, p - 5, pr, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) / (2 * I->dx[pr]);
		//	WRITE_TO_A(p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2], -1);
		//	A_value = density_t(I, p - 5, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) * avarage_velocity(I, p - 5, pr, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]) / (2 * I->dx[pr]);
		//	WRITE_TO_A(p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2], -1);
		}
	}
	return 0;
}

int DIV_density_saturation_internal_energy_avarage_velocity_backward_euler_second_separated_FDM_termogas(in *I, int p, int i, int j, int k)
{
	if (check_for_corrupt_cell(I, i, j, k)) return 1;
	double A_value, B_value;
	int pp, pr, ind_pr[3];
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
	double A_value, B_value;
	int pp, pr, ind_pr[3];
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
	for (a = 0; a < 4; a++) {
		I->B[A_IND(I, p, i, j, k)] += concentration(I, a, i, j, k) * mass_inflow_rate_func(I, a, i, j, k);
	}
	for (a = 0; a < 2; a++) {
		I->B[A_IND(I, p, i, j, k)] += saturation(I, a, i, j, k) * mass_inflow_rate_func(I, a + 4, i, j, k);
	}
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

