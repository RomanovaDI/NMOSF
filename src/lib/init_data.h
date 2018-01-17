#define AVALANCHE 0
#define TERMOGAS 1

#define DEBUG 1

typedef struct init_parameters_avalanche {
	char map_name[100];
	char region_map_name[100];
	double hight;
	int kx, ky, kz;
	int nx, ny, nz;
	int n_bl_x, n_bl_y, n_bl_z;
	int ncols;
	int nrows;
	double cellsize;
	double *mass;
	int *ind;
	int *bl_cond;
	int n_points, n_cells;
	int n_points_multipl;
	int n_cells_multipl;
	int *ind_multipl, *ind_cell_multipl;
	double interpolation, interpolation_poli;
	double *mass1;
	double nodata_value;
	int *snow_region;
	double depth;
	double density_snow, density_air;
	double pressure_atmosphere;
	double g[3];
	double k_viscosity_snow, k_viscosity_air;
	double *normal, *volume, *area;
	double dt, dx[3];
	double *B, *Aelem_csr, *B_prev;
	int *Ajptr_csr, *Aiptr_csr;
	int num_parameters;
	double shear_rate_0,limiting_viscosity_snow, flow_index, yield_stress;
	int system_dimension, system_dimension_with_boundary;
	int non_zero_elem;
	int A_ind_current;
	int flag_first_time_step;
	double end_time;
	int stencil_size;
	int *ind_boundary_cells, n_boundary_cells;
	double mass_quantity;
	int my_rank;
	int nproc;
	int *ind_proc;
	int gl_nx, gl_ny, gl_nz;
	int *gl_ind_cell_multipl;
	int x_regions, y_regions;
	int num_el_in_x_region, num_el_in_y_region;
	int max_num_el_in_x_region, max_num_el_in_y_region;
	int gl_n_cells_multipl;
	double *gl_B;
	int *ind_start_region_proc;
} in_avalanche;

typedef struct init_parameters_termogas {
	char map_name[100];
	char region_map_name[100];
	double hight;
	int kx, ky, kz;
	int nx, ny, nz;
	int n_bl_x, n_bl_y, n_bl_z;
	int ncols;
	int nrows;
	double cellsize;
	double *mass;
	int *ind;
	int *bl_cond;
	int n_points, n_cells;
	int n_points_multipl;
	int n_cells_multipl;
	int *ind_multipl, *ind_cell_multipl;
	double interpolation, interpolation_poli;
	double *mass1;
	double nodata_value;
	double *normal, *volume, *area;
	double dt, dx[3];
	double *B, *Aelem_csr, *B_prev;
	int *Ajptr_csr, *Aiptr_csr;
	int num_parameters;
	int system_dimension, system_dimension_with_boundary;
	int non_zero_elem;
	int A_ind_current;
	int flag_first_time_step;
	double end_time;
	int stencil_size;
	int *ind_boundary_cells, n_boundary_cells;
	double mass_quantity;
	double porousness;
	double density_0[2];
	double pressure_0;
	double density_coef_beta[2];
	double temperature_0;
	double density_coef_a[2];
	double R;
	double molar_weight[4];
	double permeability;
	double residual_saturation[3];
	double viscosity_coef_A[2];
	double viscosity_coef_A_gas[4];
	double viscosity_coef_B[2];
	double viscosity_coef_C_gas[4];
	double temperature_0_gas[4];
	double capillary_pressure_at_maximum_saturation[3];
	double residual_saturation_two_phase[3];
	double capillary_pressure_coef;
	double num_carbon_atoms;
	double num_hydrogen_atoms;
	double stoichiometric_coef[3];
	double threshold_temperature;
	double frequency_factor;
	double activation_temperature;
	double specific_heat[7];
	double thermal_conductivity_coef[4];
	double heat_transfer_coef;
	double tempetarure_for_calculation_internal_energy;
	double initial_enthalpy[7];
	double density_environment;
	double initial_pressure;
	double initial_temperature;
	double injection_well_pressure;
	double production_well_pressure;
	int equation_num;
	double epsilon;
	double adiabatic_exponent[6];
	double time_step;
	double volume_producted_oil_kg;
	double volume_producted_oil_m;
	int my_rank;
	int nproc;
	int *ind_proc;
	int gl_nx, gl_ny, gl_nz;
	int *gl_ind_cell_multipl;
	int x_regions, y_regions;
	int num_el_in_x_region, num_el_in_y_region;
	int max_num_el_in_x_region, max_num_el_in_y_region;
	int gl_n_cells_multipl;
	double *gl_B;
	int *ind_start_region_proc;
} in_termogas;

//typedef in_avalanche in;
typedef in_termogas in;

int set_parameters_avalanche(in *I);
int set_parameters_termogas(in *I);
