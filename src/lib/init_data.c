#include "init_data.h"
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
	strcpy(I->map_name, "maps/map_termogas.asc");
	strcpy(I->region_map_name, "maps/map_termogas_region.asc");
	//strcpy(I->map_name, "maps/map_cavity.asc");
	I->hight = 30;
	I->kx = 1;
	I->ky = 1;
	I->kz = 1;
	I->end_time = 1000;
	I->stencil_size = 2;
	I->num_parameters = 10; // 3 phase saturation, 4 components of gas concentration, pressure, temperature of porous medium, temperature of mixed flow
	I->mass_quantity = 0;
	I->dt = 1;//we need to set dt!!!
	I->porousness = 0.3;
	I->density_0[0] = 998; // water
	I->density_0[1] = 850; // oil
	I->pressure_0 = 1000000;
	I->density_coef_beta[0] = 0.0011; // water
	I->density_coef_beta[1] = 0.0008701; // oil
	I->temperature_0 = 330;
	I->density_coef_a[0] = 1400; // water
	I->density_coef_a[1] = 1300; // oil
	I->R = 8.314;
	I->molar_weight[0] = 0.02801; //N2
	I->molar_weight[1] = 0.032; //O2
	I->molar_weight[2] = 0.04401; //CO2
	I->molar_weight[3] = 0.01802; //H2O(g)
	I->permeability = 0.000000000001;
	I->residual_saturation[0] = 0.15; // water
	I->residual_saturation[1] = 0.2; // oil
	I->residual_saturation[2] = 0; // gas
	I->viscosity_coef_A[0] = 0.00000001556; // water
	I->viscosity_coef_A[1] = 0.000000383; // oil
	I->viscosity_coef_A_gas[0] = 0.00001781; // N2
	I->viscosity_coef_A_gas[1] = 0.00002018; // O2
	I->viscosity_coef_A_gas[2] = 0.0000148; // CO2
	I->viscosity_coef_A_gas[3] = 0.000018; // H2O(g)
	I->viscosity_coef_B[0] = 0.000984; // water
	I->viscosity_coef_B[1] = 0.00117; // oil
	I->viscosity_coef_C_gas[0] = 111; // N2
	I->viscosity_coef_C_gas[0] = 127; // O2
	I->viscosity_coef_C_gas[0] = 240; // CO2
	I->viscosity_coef_C_gas[0] = 0; // H2O(g)
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
	I->threshold_temperature = 400;
	I->frequency_factor = 1;
	I->activation_temperature = 400;
	I->specific_heat[0] = 4180.6; // water
	I->specific_heat[1] = 1800; // oil
	I->specific_heat[2] = 1040; // N2
	I->specific_heat[3] = 918; // O2
	I->specific_heat[4] = 820; // CO2
	I->specific_heat[5] = 2078.4; // H2O(g)
	I->specific_heat[6] = 1000; // environment
	I->thermal_conductivity_coef[0] = 0.55; // water
	I->thermal_conductivity_coef[1] = 0.128; // oil
	I->thermal_conductivity_coef[2] = 0.05; // gas
	I->thermal_conductivity_coef[3] = 2; // environment
	I->heat_transfer_coef = 0.3;
	I->tempetarure_for_calculation_internal_energy = 600;//298;
	I->initial_enthalpy[0] = -285000; // water
	I->initial_enthalpy[1] = 45000; // oil
	I->initial_enthalpy[2] = 0; // N2
	I->initial_enthalpy[3] = 0; // O2
	I->initial_enthalpy[4] = -393000; // CO2
	I->initial_enthalpy[5] = -241000; // H2O(g)
	I->initial_enthalpy[6] = 0; // environment
	I->density_environment = 2000;
	I->initial_pressure = 20000000;
	I->initial_temperature = 600;
	I->injection_well_pressure = 30000000;
	I->production_well_pressure = 100000;
	I->epsilon = 1e-10;
	I->adiabatic_exponent[0] = 1; // water
	I->adiabatic_exponent[1] = 1; // oil
	I->adiabatic_exponent[2] = 1.4; // N2
	I->adiabatic_exponent[3] = 1.4; // O2
	I->adiabatic_exponent[4] = 1.35; // CO2
	I->adiabatic_exponent[5] = 1.3; // H2O(g)
	I->time_step = 0;
	I->volume_producted_oil = 0;
#endif
	return 0;
}
