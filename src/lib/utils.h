double density(in *I, int i, int j, int k);
double velocity(in *I, int p, int i, int j, int k);
double phase_fraction(in *I, int i, int j, int k);
double pressure(in *I, int i, int j, int k);
double velocity_on_face(in *I, int p, int i, int j, int k, int s);
double density_on_face(in *I, int i, int j, int k, int s);
double strain_rate_on_face(in *I, int i, int j, int k, int s, int m, int n);
double shear_rate_on_face(in *I, int i, int j, int k, int s);
double phase_fraction_on_face(in *I, int i, int j, int k, int s);
double effective_viscosity_on_face(in *I, int i, int j, int k, int s);
double pressure_on_face(in *I, int i, int j, int k, int s);
int A_IND_S(in *I, int p, int i, int j, int k, int s);
int A_IND_S_SWITCH(in *I, int i, int j, int k, int s);
int A_IND(in *I, int p, int i, int j, int k);
int B_IND(in *I, int p, int i, int j, int k);
int AREA_IND(in *I, int i, int j, int k, int s);
int NORMAL_IND(in *I, int p, int i, int j, int k, int s);
int VOLUME_IND(in *I, int i, int j, int k);
int check_for_corrupt_cell(in *I, int i, int j, int k);
int write_to_A_csr(in *I, int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value);
int check_conservation_of_mass(in *I);