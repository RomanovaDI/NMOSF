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
	strcpy(I->map_name, "maps/termogase.asc");
	strcpy(I->region_map_name, "maps/termogase.asc");
	I->hight = 1;
	I->kx = 1;
	I->ky = 1;
	I->kz = 1;
	I->end_time = 10;
	I->stencil_size = 2;
	I->num_parameters = 10; // 3 phase saturation, 4 components of gas concentration, pressure, temperature of porous medium, temperature of mixed flow
	I->mass_quantity = 0;
	I->dt = 0.01;//we need to set dt!!!
	I->porousness = 0.3;
	I->density_0[0] = 998;
	I->density_0[1] = 850;
	I->pressure_0 = 1000000;
	I->density_coef_beta[0] = 0.0011;
	I->density_coef_beta[1] = 8701;
	I->temperature_0 = 330;
	I->density_coef_a[0] = 1400;
	I->density_coef_a[1] = 1300;
	I->R = 8.314;
	I->molar_weight[0] = 0.02801; //N2
	I->molar_weight[1] = 0.032; //O2
	I->molar_weight[2] = 0.04401; //CO2
	I->molar_weight[3] = 0.01802; //H2O
	I->permeability = 0.000000000001;
	I->relative_saturation[0] = 0.15;
	I->relative_saturation[1] = 0.2;
	I->relative_saturation[2] = 0;
	I->viscosity_coef_A[0] = 0.00000001556;
	I->viscosity_coef_A[1] = 0.000000383;
	I->viscosity_coef_A_gas[0] = 0.000018;
	I->viscosity_coef_A_gas[1] = 0.000021;
	I->viscosity_coef_A_gas[2] = 0.000018;
	I->viscosity_coef_A_gas[3] = 0.000018;
	I->viscosity_coef_B[0] = 0.000984;
	I->viscosity_coef_B[1] = 0.00117;
	I->viscosity_coef_C_gas[0] = 464;
	I->viscosity_coef_C_gas[0] = 613;
	I->viscosity_coef_C_gas[0] = -1;
	I->viscosity_coef_C_gas[0] = -89;
	I->capillary_pressure_at_maximum_saturation[0] = 55000;
	I->capillary_pressure_at_maximum_saturation[1] = 0;
	I->capillary_pressure_at_maximum_saturation[2] = 1;
	I->residual_saturation[0] = 0.3;
	I->residual_saturation[1] = 0;
	I->residual_saturation[2] = 0.3;
	I->capillary_pressure_coef = 0.5;
	I->num_carbon_atoms = 8;
	I->num_hydrogen_atoms = 18;
	I->stoichiometric_coef[0] = I->num_carbon_atoms + 0.25 * I->num_hydrogen_atoms;
	I->stoichiometric_coef[1] = I->num_carbon_atoms;
	I->stoichiometric_coef[2] = 0.25 * I->num_hydrogen_atoms;
	I->threshold_temperature = 400;
	I->frequency_factor = 1;
	I->activation_temperature = 400;
#endif
	return 0;
}
