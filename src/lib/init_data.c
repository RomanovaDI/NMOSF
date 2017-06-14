#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init_data.h"

int set_parameters(in *I)
{
	//printf("!%s!\n", I.map_name);
	strcpy(I->map_name, "maps/map_for_verification_Ab_21.asc");
	strcpy(I->region_map_name, "maps/map_for_verification_Ab_region_22.asc");
	//printf("!%s!\n", I.map_name);
	//I->region_map_name = "/home/daria/cases/NMOSF/maps/map_for_verification_Ab_region_2.asc";
	I->hight = 15;
	I->depth = 7;
	I->kx = 1;
	I->ky = 1;
	I->kz = 1;
	I->density_snow = 200;
	I->density_air = 1;
	I->pressure_atmosphere = 101325;
	I->k_viscosity_air = 0.0000148;
	I->k_viscosity_snow = 1000;
	I->flow_index = 1;
	I->yield_stress = 2000;
	I->shear_rate_0 = 0.00001;
	I->limiting_viscosity_snow = I->k_viscosity_snow * pow(I->shear_rate_0, I->flow_index - 1) + I->yield_stress / I->shear_rate_0;
	I->end_time = 0.3;
	I->g[0] = I->g[1] = I->g[2] = 0;
	I->stencil_size = 2;
	I->num_parameters = 5; // 5 = 3 components of velocity + 1 phase fraction + 1 pressure
	I->mass_quantity = 0;
	I->dt = 0.1;//we need to set dt!!!
	return 0;
}
