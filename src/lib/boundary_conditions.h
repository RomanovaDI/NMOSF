#if AVALANCHE
int SET_boundary_CONDITION_velocity_zero_gradient_on_up_and_sides_no_slip_on_low(in *I);
int SET_boundary_CONDITION_velocity_zero_gradient_on_y_sides_no_slip_on_other_upper_wall_is_mooving(in *I);
int SET_boundary_CONDITION_velocity_zero_gradient_on_y_and_x_sides_no_slip_on_other_upper_wall_is_mooving(in *I);
int SET_boundary_CONDITION_velocity_zero_gradient_on_y_and_x_and_upper_sides_no_slip_on_low(in *I);
int SET_boundary_CONDITION_velocity_fixed_value_on_all(in *I);
int SET_boundary_CONDITION_phase_fraction_zero_gradient_on_all(in *I);
int SET_boundary_CONDITION_phase_fraction_fixed_value_on_all(in *I);
int SET_boundary_CONDITION_pressure_fixed_value_on_up_and_sides_zero_gradient_on_low(in *I);
int SET_boundary_CONDITION_pressure_zero_gradient_on_all(in *I);
#endif
#if TERMOGAS
int set_injection_well(in *I, int i, int j, int k);
int set_production_well(in *I, int i, int j, int k);
int SET_boundary_CONDITION_termogas_no_boundaries_4_in_1_out(in *I);
#endif
