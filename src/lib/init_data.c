#include "init_data.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int set_parameters_avalanche(in *I)
{
#if AVALANCHE
	//strcpy(I->map_name, "maps/map_for_verification_Ab_21.asc");
	//strcpy(I->region_map_name, "maps/map_for_verification_Ab_region_22.asc");
	strcpy(I->map_name, "maps/map_cavity.asc");
	strcpy(I->region_map_name, "maps/map_cavity_region.asc");
	I->hight = 1;
	I->depth = 7;
	I->kx = 1;
	I->ky = 1;
	I->kz = 1;
	//I->density_snow = 200;
	I->density_snow = 1260;//1000;
	I->density_air = 1;
	I->pressure_atmosphere = 101325;
	I->k_viscosity_air = 0.0000148;
	I->k_viscosity_snow = 1.49;//1.49;//1000;
	I->flow_index = 1;
	I->yield_stress = 2000;
	I->shear_rate_0 = 0.00001;
	I->limiting_viscosity_snow = I->k_viscosity_snow * pow(I->shear_rate_0, I->flow_index - 1) + I->yield_stress / I->shear_rate_0;
	I->end_time = 10;
	double g = 9.81;
	double alpha = (90 * 3.14) / 180;
	I->g[0] = 0;//sin(alpha) * g;
	I->g[1] = 0;
	I->g[2] = 0;//- cos(alpha) * g;
	I->stencil_size = 2;
	I->num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	I->mass_quantity = 0;
	I->dt = 0.1;//we need to set dt!!!
#endif
	return 0;
}

int set_parameters_termogas(in *I)
{
#if TERMOGAS
	//strcpy(I->map_name, "maps/map_termogas_test_big.asc");
	strcpy(I->map_name, "maps/map_termogas_big.asc");
	//strcpy(I->map_name, "maps/map_termogas_small.asc");
	//strcpy(I->map_name, "maps/map_termogas_1d.asc");
	//strcpy(I->map_name, "maps/map_termogas_pattern.asc");
	//strcpy(I->region_map_name, "maps/map_termogas_region.asc");
	//strcpy(I->map_name, "maps/map_cavity.asc");
	I->hight = 5;//30;
	I->kx = 1;
	I->ky = 1;
	I->kz = 1;
	I->epsilon = 1e-5;
	I->end_time = 365 * 86400;//1000000;
	I->stencil_size = 2;
	I->num_parameters = 10; // 3 phase saturation, 4 components of gas concentration, pressure, temperature of porous medium, temperature of mixed flow
	I->mass_quantity = 0;
	I->dt = 1;//100.0;//10;//we need to set dt!!!
	I->porousness_value = 0.2;
	I->initial_pressure = 20000000;
	I->initial_temperature = 400;//298;
	I->density_0[0] = 998; // water
	I->density_0[1] = 850; // oil
	I->pressure_0 = 1000000;
	I->density_coef_beta[0] = 0.0011; // water
	I->density_coef_beta[1] = 0.0008701; // oil
	I->temperature_0 = 400;//330;
	I->density_coef_a[0] = 1400; // water
	I->density_coef_a[1] = 1300; // oil
	I->R = 8.314;
	I->molar_weight[0] = 0.02801; //N2
	I->molar_weight[1] = 0.032; //O2
	I->molar_weight[2] = 0.04401; //CO2
	I->molar_weight[3] = 0.01802; //H2O(g)
	I->molar_weight[4] = 0.11423; //C8H18
	I->permeability = 1e-14;
	I->residual_saturation[0] = 0.15; // water
	I->residual_saturation[1] = 0.2; // oil
	I->residual_saturation[2] = 0.05;//I->epsilon; // gas
	I->viscosity_coef_A[0] = 0.00000001556; // water
	I->viscosity_coef_A[1] = 0.000000383; // oil
	I->viscosity_coef_A_gas[0] = 0.00001781; // N2
	I->viscosity_coef_A_gas[1] = 0.00002018; // O2
	I->viscosity_coef_A_gas[2] = 0.0000148; // CO2
	I->viscosity_coef_A_gas[3] = 0.000018; // H2O(g)
	I->viscosity_coef_B[0] = 0.0000984; // water
	I->viscosity_coef_B[1] = 0.000117;// - 0.00007; // oil
	I->viscosity_coef_C_gas[0] = 111; // N2
	I->viscosity_coef_C_gas[1] = 127; // O2
	I->viscosity_coef_C_gas[2] = 240; // CO2
	I->viscosity_coef_C_gas[3] = 0; // H2O(g)
	I->temperature_0_gas[0] = 300.55; // N2
	I->temperature_0_gas[1] = 292.25; // O2
	I->temperature_0_gas[2] = 293.15; // CO2
	I->temperature_0_gas[3] = 330; // H2O(g)
	I->capillary_pressure_at_maximum_saturation[0] = 55000; // water
	I->capillary_pressure_at_maximum_saturation[1] = 0; // oil
	I->capillary_pressure_at_maximum_saturation[2] = 1; // gas
	I->residual_saturation_two_phase[0] = 0.3; // in system "oil - water"
	I->residual_saturation_two_phase[1] = 0; // in system "oil - oil"
	I->residual_saturation_two_phase[2] = 0.3; // in system "oil - gas"
	I->capillary_pressure_coef = 0.5;
	I->num_carbon_atoms = 8;
	I->num_hydrogen_atoms = 18;
	I->stoichiometric_coef[0] = I->num_carbon_atoms + 0.25 * I->num_hydrogen_atoms;
	I->stoichiometric_coef[1] = I->num_carbon_atoms;
	I->stoichiometric_coef[2] = 0.25 * I->num_hydrogen_atoms;
	I->stoichiometric_coef_before[0] = 0;//N2
	I->stoichiometric_coef_before[1] = 12.5;//O2
	I->stoichiometric_coef_before[2] = 0;//CO2
	I->stoichiometric_coef_before[3] = 0;//H2O(g)
	I->stoichiometric_coef_before[4] = 0;//water
	I->stoichiometric_coef_before[5] = 1;//oil
	I->stoichiometric_coef_after[0] = 0;//N2
	I->stoichiometric_coef_after[1] = 0;//O2
	I->stoichiometric_coef_after[2] = 8;//CO2
	I->stoichiometric_coef_after[3] = 4.5;//H2O(g)
	I->stoichiometric_coef_after[4] = 4.5;//water
	I->stoichiometric_coef_after[5] = 0;//oil
	I->threshold_temperature = 401;
	I->frequency_factor = 1000;
	I->activation_temperature = 402;//I->initial_temperature;//400;
	I->stoichiometric_coef_activ = 1;
	I->pressure_activ = I->initial_pressure;
	I->specific_heat[0] = 4180.6; // water
	I->specific_heat[1] = 2500;//1800; // oil
	I->specific_heat[2] = 1040; // N2
	I->specific_heat[3] = 918; // O2
	I->specific_heat[4] = 820; // CO2
	I->specific_heat[5] = 2078.4; // H2O(g)
	I->specific_heat[6] = 1000; // environment
	I->thermal_conductivity_coef[0] = 0.55; // water
	I->thermal_conductivity_coef[1] = 0.128; // oil
	I->thermal_conductivity_coef[2] = 0.05; // gas
	I->thermal_conductivity_coef[3] = 2; // environment
	I->heat_transfer_coef = 1100;
	I->tempetarure_for_calculation_internal_energy = I->initial_temperature;//298;
	I->initial_enthalpy[0] = -285800; // water
	I->initial_enthalpy[1] = -249950; // oil
	I->initial_enthalpy[2] = 0; // N2
	I->initial_enthalpy[3] = 0; // O2
	I->initial_enthalpy[4] = -393500; // CO2
	I->initial_enthalpy[5] = -241800; // H2O(g)
	I->initial_enthalpy[6] = 0; // environment
	I->density_environment = 2000;
	//I->injection_well_pressure = 30000000;
	I->injection_well_pressure = 30000000;
	I->injection_well_temperature = 600;//I->initial_temperature;//600;
	I->production_well_pressure = 10000000;
	I->adiabatic_exponent[0] = 1; // water
	I->adiabatic_exponent[1] = 1; // oil
	I->adiabatic_exponent[2] = 1.4; // N2
	I->adiabatic_exponent[3] = 1.4; // O2
	I->adiabatic_exponent[4] = 1.35; // CO2
	I->adiabatic_exponent[5] = 1.3; // H2O(g)
	I->time_step = 0;
	I->volume_producted_oil_kg = 0;
	I->volume_producted_oil_m = 0;
	I->volume_producted_fluid_m = 0;
	I->volume_injected_fluid_m = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &I->my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &I->nproc);
	I->x_regions = 1;
	I->y_regions = 1;
//	I->x_regions = 1;
//	I->y_regions = 2;
	if (I->x_regions * I->y_regions != I->nproc) {
		printf("Number of processes is not equal to number of subdomains\n");
		return 1;
	}
	I->time = 0;
	I->dependent_variables = 71;
	I->nan_flag = 0;
	I->negative_num_flag = 0;
#define TIME_STEP 0
#define SECOND 1
#define MINUTE 60
#define HOUR 3600
#define DAY 86400
#define MONTH 2592000
#define YEAR 31536000
	I->units_for_write_interval = DAY; // {TIME_STEP, SECOND, MINUTE, HOUR, DAY, MONTH, YEAR}
	I->write_interval = 1;
	I->written_step = -1;
	I->courant_number = 2;
#endif
	return 0;
}
