#include "init_data.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void)
{
	in II;
	in *I = &II;
	read_I_termogas(I);
	int num_res = 43984;
	int num_plots = 5;
	for (int i = 0; i <= num_res; i += num_res / num_plots) {
		char file_name[14 + 20];
		sprintf(file_name, "result/map%d.vtk", i / (num_res / num_plots));
		FILE *f = fopen(file_name, "w");
		sprintf(file_name, "result/s_o%d.dat", i);
		FILE *f1 = fopen(file_name, "w");
		char s[100];
		int j_pointer = 0;
		int j;
		for (j = j_pointer; j < j_pointer + 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		double tmp;
		if (sscanf(s, "oil displacement, time = %lf", &tmp))
			return 1;
		fprintf(f1, "#time = %lf\n", tmp);
		tmp = j_pointer + 2 +
			I->n_points_multipl * ((int) (I->hight / (I->cellsize / (double) I->kz)) + 1) +
			I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) +
			I->n_cells * I->kx * I->ky * I->kz * (int) (I->hight / I->cellsize) +
			I->n_cells_multipl * I->nz;
		for (j = j_pointer; j < (int) tmp; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		for (j = j_pointer; j < j_pointer + 10 * (I->nx * I->ny * I->nz + 2); j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		for (j = j_pointer; j < j_pointer + (I->nx / 2) * I->ny + I->ny / 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		sscanf(s, "%lf", &tmp);
		double step = 0;
		fprintf(f1, "%lf\t%lf\n", step, tmp);
		for (int k = 0; k < I->nx / 2; k++) {
			for (j = j_pointer; j < j_pointer + I->ny + 1; j)
				if (fgets(s, 100, f) == NULL)
					return 1;
			j_pointer = j;
			step += I->cellsize * sqrt(2);
			sscanf(s, "%lf", &tmp);
			fprintf(f1, "%lf\t%lf\n", step, tmp);
		}
		fclose(f1);
		sprintf(file_name, "result/T%d.dat", i / (num_res / num_plots));
		f1 = fopen(file_name, "w");
		for (j = j_pointer; j < j_pointer + I->nx * I->ny * I->nz + 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		for (j = j_pointer; j < j_pointer + (I->nx / 2) * I->ny + I->ny / 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		sscanf(s, "%lf", &tmp);
		step = 0;
		fprintf(f1, "%lf\t%lf\n", step, tmp);
		for (int k = 0; k < I->nx / 2; k++) {
			for (j = j_pointer; j < j_pointer + I->ny + 1; j)
				if (fgets(s, 100, f) == NULL)
					return 1;
			j_pointer = j;
			step += I->cellsize * sqrt(2);
			sscanf(s, "%lf", &tmp);
			fprintf(f1, "%lf\t%lf\n", step, tmp);
		}
		fclose(f1);
		sprintf(file_name, "result/omega%d.dat", i / (num_res / num_plots));
		f1 = fopen(file_name, "w");
		for (j = j_pointer; j < j_pointer + I->nx * I->ny * I->nz + 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		for (j = j_pointer; j < j_pointer + (I->nx / 2) * I->ny + I->ny / 2; j++)
			if (fgets(s, 100, f) == NULL)
				return 1;
		j_pointer = j;
		sscanf(s, "%lf", &tmp);
		step = 0;
		fprintf(f1, "%lf\t%lf\n", step, tmp);
		for (int k = 0; k < I->nx / 2; k++) {
			for (j = j_pointer; j < j_pointer + I->ny + 1; j)
				if (fgets(s, 100, f) == NULL)
					return 1;
			j_pointer = j;
			step += I->cellsize * sqrt(2);
			sscanf(s, "%lf", &tmp);
			fprintf(f1, "%lf\t%lf\n", step, tmp);
		}
		fclose(f1);
		fclose(f);
	}
	return 0;
}
