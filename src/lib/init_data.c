#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init_data.h"

int set_parameters(in *I)
{
	strcpy(I->map_name, "maps/map_for_verification_Ab_21.asc");
	strcpy(I->region_map_name, "maps/map_for_verification_Ab_region_22.asc");
	strcpy(I->map_name, "maps/map_cavity.asc");
	strcpy(I->region_map_name, "maps/map_cavity_region.asc");
	I->hight = 0.1;
	I->depth = 7;
	I->kx = 20;
	I->ky = 1;
	I->kz = 20;
	//I->density_snow = 200;
	I->density_snow = 1000;
	I->density_air = 1;
	I->pressure_atmosphere = 101325;
	I->k_viscosity_air = 0.0000148;
	I->k_viscosity_snow = 1.49;//1000;
	I->flow_index = 1;
	I->yield_stress = 2000;
	I->shear_rate_0 = 0.00001;
	I->limiting_viscosity_snow = I->k_viscosity_snow * pow(I->shear_rate_0, I->flow_index - 1) + I->yield_stress / I->shear_rate_0;
	I->end_time = 100;
	double g = 9.81;
	double alpha = (90 * 3.14) / 180;
	I->g[0] = 0;//sin(alpha) * g;
	I->g[1] = 0;
	I->g[2] = 0;//- cos(alpha) * g;
	I->stencil_size = 2;
	I->num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	I->mass_quantity = 0;
	I->dt = 0.1;//we need to set dt!!!
	return 0;
}
