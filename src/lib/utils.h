int A_IND_S(in *I, int p, int i, int j, int k, int s);
int A_IND_S_SWITCH(in *I, int i, int j, int k, int s);
int A_IND(in *I, int p, int i, int j, int k);
int GL_A_IND(in *I, int p, int i, int j, int k);
int B_IND(in *I, int p, int i, int j, int k);
int PARAM_IND(in *I, int p, int i, int j, int k);
int POR_IND(in *I, int i, int j, int k);
int INT_NEI_IND(in *I, int p, int i, int j, int k);
int POR_INT_NEI_IND(in *I, int i, int j, int k);
int AREA_IND(in *I, int i, int j, int k, int s);
int NORMAL_IND(in *I, int p, int i, int j, int k, int s);
int VOLUME_IND(in *I, int i, int j, int k);
int check_for_corrupt_cell(in *I, int i, int j, int k);
int well(in *I, int i, int j, int k);
int production_well(in *I, int i, int j, int k);
int injection_well(in *I, int i, int j, int k);
int write_to_A_csr(in *I, int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value);
int cell_of_computation_domain(in *I, int i, int j, int k);
int internal_cell(in *I, int i, int j, int k);
int boundary_cell(in *I, int i, int j, int k);
int count_neighbor_internal_cells(in *I, int i, int j, int k);
int count_second_order_neighbor_internal_cells(in *I, int i, int j, int k);
int count_other_corner_neighbor_internal_cells(in *I, int i, int j, int k);
double max2(double a, double b);
double min2(double a, double b);
double max3(double a, double b, double c);
double min3(double a, double b, double c);
int nan_check(in *I, double value, int p, int pr, int i, int j, int k, const char func_name[100]);
int negative_or_zero_num_check(in *I, double value, int p, int pr, int i, int j, int k, const char func_name[100]);
#if AVALANCHE
int check_conservation_of_mass(in *I);
double density(in *I, int i, int j, int k);
double velocity(in *I, int p, int i, int j, int k);
double phase_fraction(in *I, int i, int j, int k);
double pressure(in *I, int i, int j, int k);
double velocity_on_face(in *I, int p, int i, int j, int k, int s);
double density_on_face(in *I, int i, int j, int k, int s);
double strain_rate_on_face(in *I, int i, int j, int k, int s, int m, int n);
double shear_rate_on_face(in *I, int i, int j, int k, int s);
double strain_rate(in *I, int i, int j, int k, int m, int n);
double shear_rate(in *I, int i, int j, int k);
double phase_fraction_on_face(in *I, int i, int j, int k, int s);
double effective_viscosity_on_face(in *I, int i, int j, int k, int s);
double pressure_on_face(in *I, int i, int j, int k, int s);
int barotropy_pressure(in *I);
int barotropy_density(in *I);
#endif
#if TERMOGAS
double saturation(in *I, int p, int i, int j, int k);
double concentration(in *I, int p, int i, int j, int k);
double pressure(in *I, int i, int j, int k);
double temperature_flow(in *I, int i, int j, int k);
double temperature_environment(in *I, int i, int j, int k);
double density_t_func(in *I, int p, int i, int j, int k);
double density_t(in *I, int p, int i, int j, int k);
double two_phase_relative_permeability_func(in *I, int p, int pr, int i, int j, int k);
double two_phase_relative_permeability(in *I, int p, int pr, int i, int j, int k);
double relative_permeability_func(in *I, int p, int i, int j, int k);
double relative_permeability(in *I, int p, int i, int j, int k);
double viscosity_gas_component_func(in *I, int p, int i, int j, int k);
double viscosity_gas_component(in *I, int p, int i, int j, int k);
double molar_fraction_func(in *I, int p, int i, int j, int k);
double molar_fraction(in *I, int p, int i, int j, int k);
double viscosity_func(in *I, int p, int i, int j, int k);
double viscosity(in *I, int p, int i, int j, int k);
double Darsi_M_coef_phases_func(in *I, int p, int i, int j, int k);
double Darsi_M_coef_phases(in *I, int p, int i, int j, int k);
double Darsi_M_coef_func(in *I, int i, int j, int k);
double Darsi_M_coef(in *I, int i, int j, int k);
double capillary_pressure_derivative_by_saturation_func(in *I, int p, int i, int j, int k);
double capillary_pressure_derivative_by_saturation(in *I, int p, int i, int j, int k);
double grad_pressure_func(in *I, int pr, int i, int j, int k);
double grad_pressure(in *I, int pr, int i, int j, int k);
double grad_saturation_func(in *I, int p, int pr, int i, int j, int k);
double grad_saturation(in *I, int p, int pr, int i, int j, int k);
double coef_grad_saturation_func(in *I, int p, int pr, int i, int j, int k);
double coef_grad_saturation(in *I, int p, int pr, int i, int j, int k);
double avarage_velocity_func(in *I, int p, int pr, int i, int j, int k);
double avarage_velocity(in *I, int p, int pr, int i, int j, int k);
double rate_of_reaction_coef_func(in *I, int i, int j, int k);
double rate_of_reaction_coef(in *I, int i, int j, int k);
double rate_of_reaction_func(in *I, int i, int j, int k);
double rate_of_reaction(in *I, int i, int j, int k);
double mass_inflow_rate_func(in *I, int p, int i, int j, int k);
double mass_inflow_rate(in *I, int p, int i, int j, int k);
double chemical_reaction_heat_flow_func(in *I, int i, int j, int k);
double chemical_reaction_heat_flow(in *I, int i, int j, int k);
double density_derivative_by_pressure(in *I, int p, int i, int j, int k);
double Darsi_A_coef_func(in *I, int i, int j, int k);
double Darsi_A_coef(in *I, int i, int j, int k);
double thermal_conductivity_func(in *I, int i, int j, int k);
double thermal_conductivity(in *I, int i, int j, int k);
double convection_speed_func(in *I, int pr, int i, int j, int  k);
double convection_speed(in *I, int pr, int i, int j, int  k);
int check_sum(in *I);
int check_convection_speed(in *I);
int check_values_in_B(in *I);
int check_values_in_B_prev(in *I);
int check_for_symmetry(in *I);
int print_oil_production(in *I);
double avarage_velocity_global(in *I);
int save_I_termogas(in *I);
int read_I_termogas(in *I);
int calculate_disparity_pressure(in *I);
#endif
