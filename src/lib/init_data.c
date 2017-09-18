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
	I->density_water = 998;
	I->density_oil = 850;
	I->pressure_0 = 1000000;
	I->beta_water = 0.0011;
	I->beta_oil = 8701;
	I->temperature_0 = 330;
	I->a_water = 1400;
	I->a_oil = 1300;
	I->R = 8.314;
	I->m_weight[0] = 0.02801; //N2
	I->m_weight[1] = 0.032; //O2
	I->m_weight[2] = 0.04401; //CO2
	I->m_weight[3] = 0.01802; //H2O
	I->permeability = 0.000000000001;
	I->relative_saturation[0] = 0.15;
	I->relative_saturation[1] = 0.2;
	I->relative_saturation[2] = 0;
#endif
	return 0;
}
