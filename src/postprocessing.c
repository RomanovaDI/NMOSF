#include "init_data.h"
#include "read_files.h"
#include "mesh_operations.h"
#include "utils.h"
#include "vtk_map_functions.h"
#include "boundary_conditions.h"
#include "initial_conditions.h"
#include "array_functions.h"
#if AVALANCHE
#include "x_crank_nikolson_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_VOF_avalanche.h"
#include "x_forward_euler_second_combined_FDM_avalanche.h"
#include "x_backward_euler_second_combined_VOF_avalanche.h"
#include "x_backward_euler_second_combined_FDM_avalanche.h"
#include "t_second_combined_VOF_avalanche.h"
#endif
#if TERMOGAS
#include "t_second_combined_FDM_termogas.h"
#include "x_backward_euler_second_combined_FDM_termogas.h"
#include "t_second_separated_FDM_termogas.h"
#include "x_backward_euler_second_separated_FDM_termogas.h"
#endif
#include "t_test.h"
#include "create_matrix.h"
#include "matrix_functions.h"
#include "slu_ddefs.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#define SET_CONDITION(type, object, mode) SET_##type##_CONDITION_##object##_##mode(I)

int read_map(int map_num, in *I)
{
	char file_name[14 + 20];
	sprintf(file_name, "result/map%d.vtk", map_num);
	FILE *f = fopen(file_name, "r");
	char s[100], s1[100];
	for (int j = 0; j < 2; j++) {
		if (fgets(s, 100, f) == NULL)
			return 1;
	}
	double tmp;
	sscanf(s, "%s %lf", s1, &I->time);
	int tmp_j = 2 +
		I->n_points_multipl * ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1) + 1 +
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) + 1 +
		I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) + 1 + 1;
	for (int j = 0; j < tmp_j; j++)
		if (fgets(s, 100, f) == NULL)
			return 1;
	int p = 0;
	while (p < I->num_parameters) {
		if (fgets(s, 100, f) == NULL)
			return 1;
		char s2[100], s3[100];
		sscanf(s, "%s %s %s", s1, s2, s3);
		char ss[10][100] = {
			"concentration_0",
			"concentration_1",
			"concentration_2",
			"concentration_3",
			"pressure",
			"saturation_0",
			"saturation_1",
			"saturation_2",
			"temperature_flow",
			"temperature_environment"};
		int flag = 0;
		for (int i = 0; i < I->num_parameters; i++)
			if (! strcmp(s2, ss[i]))
				flag = 1;
		if (flag) {
			if (fgets(s, 100, f) == NULL)
				return 1;
			for (int k = 0; k < I->gl_nz; k++) {
				for (int i = 0; i < I->gl_nx; i++) {
					for (int j = 0; j < I->gl_ny; j++) {
						if (fgets(s, 100, f) == NULL)
							return 1;
						sscanf(s, "%lf", &I->B[A_IND(I, p, i, j, k)]);
						if ((p == 8) || (p == 9))
							I->B[A_IND(I, p, i, j, k)] -= I->initial_temperature;
					}
				}
			}
			p++;
		} else {
			for (int i = 0; i < I->gl_n_cells_multipl * I->gl_nz + 1; i++)
				if (fgets(s, 100, f) == NULL)
					return 1;
		}
	}
	fclose(f);
	return 0;
}

int print_injection_production_param(int p, int step_num, in *I)
{
	char file_name[14 + 20];
	sprintf(file_name, "result/parameter%d_step%d.dat", p, step_num);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "#time = %lf\n", I->time);
	double step = 0;
	for (int i = I->gl_nx / 2; i < I->gl_nx; i++) {
		fprintf(f, "%20.20lf %20.20lf\n", step, I->B[A_IND(I, p, i, i, 0)]);
		step += I->cellsize * sqrt(2);
	}
	fclose(f);
	return 0;
}

int print_injection_production_rate_of_reaction(int step_num, in *I)
{
	char file_name[100];
	sprintf(file_name, "result/rate_of_reaction_step%d.dat", step_num);
	FILE *f = fopen(file_name, "w");
	fprintf(f, "#time = %lf\n", I->time);
	double step = 0;
	for (int i = I->gl_nx / 2; i < I->gl_nx; i ++) {
		fprintf(f, "%20.20lf %20.20lf\n", step, rate_of_reaction(I, i, i, 0));
		step += I->cellsize * sqrt(2);
	}
	fclose(f);
	return 0;
}

int print_symmetrical_case(int step_num, in *I)
{
	int k = 0;
	for (int i = 0; i < I->gl_nx / 2; i++) {
		for (int j = 0; j < I->gl_ny / 2; j++) {
			for (int p = 0; p < I->num_parameters; p++) {
				I->B_prev[B_IND(I, p, i, I->gl_ny - 1 - j, k)] = I->B_prev[B_IND(I, p, i, j, k)];
				I->B_prev[B_IND(I, p, I->gl_nx - 1 - i, j, k)] = I->B_prev[B_IND(I, p, i, j, k)];
				I->B_prev[B_IND(I, p, I->gl_nx - 1 - i, I->gl_ny - 1 - j, k)] = I->B_prev[B_IND(I, p, i, j, k)];
			}
		}
	}
	I->time_step = -step_num;
	if (print_vtk(I))
		return 1;
	return 0;
}

int main(int argc, char **argv)
{
	in II;
	in *I = &II;
	//read_I_termogas(I);
	MPI_Init(&argc, &argv);
	I->flag_first_time_step = 1;
	if (set_parameters_termogas(I)) return 1;
	MPI_Barrier(MPI_COMM_WORLD);
	if (read_asc_and_declare_variables(I)) return 1;
	if (I->my_rank == 0)
		if (do_interpolation(I)) return 1;
	if (I->nproc > 1)
		if (share_dx(I)) return 1;
	if (do_decomposition(I)) return 1;
	if (make_boundary(I)) return 1;
	I->system_dimension = I->n_cells_multipl * I->nz * I->num_parameters;
	I->system_dimension_with_boundary = I->n_boundary_cells * (I->nz + 2 * I->stencil_size) * I->num_parameters;
	if ((I->B = (double *) malloc(I->system_dimension * sizeof(double))) == NULL) {
		printf("Memory error in function %s.\n", __func__);
		return 1;
	}
	if ((I->B_prev = (double *) malloc(I->system_dimension_with_boundary * sizeof(double))) == NULL) {
		printf("Memory error in function %s.\n", __func__);
		return 1;
	}
	int num_res = 268910;
	int num_plots = 5;
	for (int i = 0; i <= num_res; i += num_res / num_plots) {
		printf("%d\n", i);
		read_map(i, I);
		if (set_array_of_parameters_termogas(I)) return 1;
		I->flag_first_time_step = 0;
		if (write_B_to_B_prev(I)) return 1;
		SET_CONDITION(boundary, termogas, no_boundaries_4_in_1_out);
		for (int j = 0; j < I->num_parameters; j++)
			print_injection_production_param(j, i / (num_res / num_plots), I);
		print_injection_production_rate_of_reaction(i / (num_res / num_plots), I);
		if (print_symmetrical_case(i, I))
			return 1;
	}
	free(I->B);
	free(I->B_prev);
	free(I->array_of_parameters);
	MPI_Finalize();
	return 0;
}
