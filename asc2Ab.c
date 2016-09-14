#include<stdio.h>
#include<stdlib.h>
#include<string.h>

//asc2vtk.out [map.asc] [hight] [x koef] [y koef] [z koef] [region_map.asc]


int read_asc_and_declare_variables(int argc, char **argv);
int do_interpolation(int argc, char **argv);
int print_vtk(int argc, char **argv);
int free_massives(int argc, char **argv);

float hight;
int ncols;
int nrows;
float cellsize;
float *mass;
int *ind;
int *bl_cond;
int n_points, n_cells;
int n_points_multipl;
int *ind_multipl;
int count_ind_multipl;
float interpolation, interpolation_poli;
float *mass1;
float nodata_value;

int read_asc_and_declare_variables(int argc, char **argv)
{
	extern float hight;
	extern int ncols;
	extern int nrows;
	extern float cellsize;
	extern float *mass;
	extern int *ind;
	extern int *bl_cond;
	extern int n_points, n_cells;
	extern int n_points_multipl;
	extern int *ind_multipl;
	extern int count_ind_multipl;
	extern float interpolation, interpolation_poli;
	extern float *mass1;
	extern float nodata_value;
	
	FILE *f = fopen(argv[1],"r");
	if (f == NULL) {
		printf("No such file\n");
		return 1;;
	}

	FILE *f1 = fopen("map.txt", "w");
	int i;
	while ((i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);

	f = fopen(argv[6],"r");
	if (f == NULL) {
		printf("No such file\n");
		goto exit_without_massives;
	}

	f1 = fopen("regions_map.txt", "w");
	while ((i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);

	f = fopen("map.txt", "r");
	f1 = fopen("regions_map.txt", "r");
	int err;
	char str[20];

	hight = atof(argv[2]);
	
	if ((err = fscanf(f, "%s %d", str, &ncols)) == EOF) {
		printf("Error file of map\n");
		return 1;
	}
	if (strcmp(str, "ncols") != 0) return 1;

	if ((err = fscanf(f, "%s %d", str, &nrows)) == EOF) return 1;
	if (strcmp(str, "nrows") != 0) return 1;

	float xllcorner;
	if ((err = fscanf(f, "%s %f", str, &xllcorner)) == EOF) return 1;
	if (strcmp(str, "xllcorner") != 0) return 1;
	
	float yllcorner;
	if ((err = fscanf(f, "%s %f", str, &yllcorner)) == EOF) return 1;
	if (strcmp(str, "yllcorner") != 0) return 1;
	
	if ((err = fscanf(f, "%s %f", str, &cellsize)) == EOF) return 1;
	if (strcmp(str, "cellsize") != 0) return 1;

	if ((err = fscanf(f, "%s %f", str, &nodata_value)) == EOF) return 1;
	if (strcmp(str, "NODATA_value") != 0) return 1;

	int ncols1;
	if ((err = fscanf(f1, "%s %d", str, &ncols1)) == EOF) return 1;;
	if (strcmp(str, "ncols") != 0) goto err_file;

	int nrows1;
	if ((err = fscanf(f1, "%s %d", str, &nrows1)) == EOF) return 1;;
	if (strcmp(str, "nrows") != 0) goto err_file;

	float xllcorner1;
	if ((err = fscanf(f1, "%s %f", str, &xllcorner1)) == EOF) return 1;
	if (strcmp(str, "xllcorner") != 0) goto err_file;
	
	float yllcorner1;
	if ((err = fscanf(f1, "%s %f", str, &yllcorner1)) == EOF) return 1;
	if (strcmp(str, "yllcorner") != 0) goto err_file;
	
	float cellsize1;
	if ((err = fscanf(f1, "%s %f", str, &cellsize1)) == EOF) return 1;
	if (strcmp(str, "cellsize") != 0) goto err_file;

	float nodata_value1;
	if ((err = fscanf(f1, "%s %f", str, &nodata_value1)) == EOF) return 1;
	if (strcmp(str, "NODATA_value") != 0) goto err_file;

	if (cellsize != cellsize1) {
		printf("Cellsize in both maps need to be the same\n");
		return 1;
	}

	if ((cellsize - (int) cellsize != 0) || (cellsize1 - (int) cellsize1 != 0)) {
		printf("In this vercion the value cellsize need to be integer\n");
		return 1;
	}

	if ((((int) fabs(xllcorner - xllcorner1)) % (int) cellsize) || ((int) fabs(yllcorner - yllcorner1) % (int) cellsize) ||
		(fabs(xllcorner - xllcorner1) - (int) fabs(xllcorner - xllcorner1) != 0) || (fabs(yllcorner - yllcorner1) - (int) fabs(yllcorner - yllcorner1) != 0)) {
		printf("Difference between xllcorners of maps and yllcorners need to aligned to cellsize\n");
		return 1;
	}

	if ((mass = (float *) malloc(ncols * nrows * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((ind = (int *) malloc(ncols * nrows * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	if ((bl_cond = (int *) malloc((ncols - 1) * (nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	for (i = 0; i < ncols * nrows; i++) {
		if (err = fscanf(f, "%f", &mass[i]) == EOF)
			return 1;
		ind[i] = -1;
		if (i < (ncols - 1) * (nrows - 1))
			bl_cond[i] = -1;
	}
	fclose(f);
	remove("map.txt");

	int j;
	n_points = 0;
	n_cells = 0;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if ((i - 1 >= 0) && (j + 1 < ncols) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[(i - 1) * ncols + j] != nodata_value) &&
				(mass[i * ncols + j + 1] != nodata_value) &&
				(mass[(i - 1) * ncols + j + 1] != nodata_value)) {
					ind[i * ncols +j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((j + 1 < ncols) && (i + 1 < nrows) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[i * ncols + j + 1] != nodata_value) &&
				(mass[(i + 1) * ncols + j] != nodata_value) &&
				(mass[(i + 1) * ncols + j + 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((i + 1 < nrows) && (j - 1 >= 0) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[(i + 1) * ncols + j] != nodata_value) &&
				(mass[i * ncols + j - 1] != nodata_value) &&
				(mass[(i + 1) * ncols + j - 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(mass[i * ncols + j] != nodata_value) &&
				(mass[i * ncols + j - 1] != nodata_value) &&
				(mass[(i - 1) * ncols + j] != nodata_value) &&
				(mass[(i - 1) * ncols + j - 1] != nodata_value)) {
					ind[i * ncols + j] = n_points;
					n_points++;
					goto end_points;
			}
			end_points:
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(ind[i * ncols + j] != -1) &&
				(ind[i * ncols + j - 1] != -1) &&
				(ind[(i - 1) * ncols + j - 1] != -1) &&
				(ind[(i - 1) * ncols + j] != -1)) {
					bl_cond[(i - 1) * (ncols - 1) + j - 1] = 1;
					n_cells++;
			}
		}
	}
	return 0;
}

int do_interpolation(int argc, char **argv)
{
	extern float hight;
	extern int ncols;
	extern int nrows;
	extern float cellsize;
	extern float *mass;
	extern int *ind;
	extern int *bl_cond;
	extern int n_points, n_cells;
	extern int n_points_multipl;
	extern int *ind_multipl;
	extern int count_ind_multipl;
	extern float interpolation, interpolation_poli;
	extern float *mass1;
	extern float nodata_value;

	int i, j, k, l, m, n, o, p, q, a = 0;
	n_points_multipl = 0;
	for (i = 0; i < nrows - 1; i++) {
		for (j = 0; j < ncols - 1; j++) {
			if (bl_cond[i * (ncols - 1) + j] == 1) {
				if ((j - 1 < 0) || (bl_cond[i * (ncols - 1) + j - 1] != 1)) {
					n_points_multipl += atoi(argv[4]) + 1;
					if ((i - 1 >= 0) && (bl_cond[(i - 1) * (ncols - 1) + j] == 1)) {
						n_points_multipl -= 1;
					}
				}
				if ((i - 1 < 0) || (bl_cond[(i - 1) * (ncols - 1) + j] != 1)) {
					n_points_multipl += atoi(argv[3]);
					if ((i - 1 >= 0) && (j + 1 < ncols - 1) && (bl_cond[(i - 1) * (ncols - 1) + j + 1] == 1)) {
						n_points_multipl -= 1;
					}
				}
				n_points_multipl += atoi(argv[3]) * atoi(argv[4]);
			}
		}
	}

	if ((ind_multipl = (int *) malloc(ncols * atoi(argv[3]) * nrows * atoi(argv[4]) * sizeof(int))) == NULL) {
		printf("Memory error\n");
		return 1;
	}

	count_ind_multipl = 0;
	if ((mass1 = (float *) malloc(ncols * atoi(argv[3]) * nrows * atoi(argv[4]) * sizeof(float))) == NULL) {
		printf("Memory error\n");
		return 1;
	}
	for (i = 0; i < nrows; i++) {
		for (k = 0; k < atoi(argv[4]); k++) { 
			for (j = 0; j < ncols; j++) {
				for (l = 0; l < atoi(argv[3]); l++) {
					mass1[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = nodata_value;
					ind_multipl[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = -1;
					if (ind[i * ncols + j] != -1) {
						if ((k == 0) && (l == 0)) {
							mass1[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = mass[i * ncols + j];
							ind_multipl[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = count_ind_multipl;
							count_ind_multipl++;
						} else {
							if (((k == 0) && (j + 1 < ncols) && (ind[i * ncols + j + 1] != -1)) ||
								((l == 0) && (i + 1 < nrows) && (ind[(i + 1) * ncols + j] != -1)) ||
								((i + 1 < nrows) && (j + 1 < ncols) && (ind[(i + 1) * ncols + j] != -1) &&
								 (ind[i * ncols + j + 1] != -1) && (ind[(i + 1) * ncols + j + 1] != -1))) {
									interpolation = 0;
									for (n = 0; n < nrows; n++) {
										for (m = 0; m < ncols; m++) {
											if (ind[n * ncols + m] != -1) {
												interpolation_poli = mass[n * ncols + m];
												for (o = 0; o < nrows; o++) {
													if (o != n) {
														interpolation_poli *= (i + k / atof(argv[4]) - o) / (n - o);
													}
												}
												for (p = 0; p < ncols; p ++) {
													if (p != m) {
														interpolation_poli *= (j + l / atof(argv[3]) - p) / (m - p);
													}
												}
												interpolation += interpolation_poli;
											}
										}
									}
									mass1[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = interpolation;
									ind_multipl[i * atoi(argv[4]) * ncols * atoi(argv[3]) + k * ncols * atoi(argv[3]) + j * atoi(argv[3]) + l] = count_ind_multipl;
									count_ind_multipl++;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int print_vtk(int argc, char **argv)
{
	extern float hight;
	extern int ncols;
	extern int nrows;
	extern float cellsize;
	extern float *mass;
	extern int *ind;
	extern int *bl_cond;
	extern int n_points, n_cells;
	extern int n_points_multipl;
	extern int *ind_multipl;
	extern int count_ind_multipl;
	extern float interpolation, interpolation_poli;
	extern float *mass1;

	FILE *f = fopen("map.vtk", "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "slope\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

	int i, j, k, a;
	fprintf(f, "POINTS %d float\n", n_points_multipl * ((int) (hight / cellsize) * atoi(argv[5]) + 1));
	for (k = 0; k <= (int) (hight / cellsize) * atoi(argv[5]); k++) {
		for (i = 0; i < nrows * atoi(argv[4]); i ++) {
			for (j = 0; j < ncols * atoi(argv[3]); j++) {
				if (ind_multipl[i * ncols * atoi(argv[3]) + j] != -1)
					fprintf(f, "%f %f %f\n", j * cellsize / atof(argv[3]), i * cellsize / atof(argv[4]),
						mass1[i * ncols * atoi(argv[3]) + j] + k * cellsize / atof(argv[5]));
			}
		}
	}
	
	a = 0;
	fprintf(f, "CELLS %d %d\n", n_cells * atoi(argv[3]) * atoi(argv[4]) * atoi(argv[5]) * (int) (hight / cellsize) * 5,
		n_cells * atoi(argv[3]) * atoi(argv[4]) * atoi(argv[5]) * (int) (hight / cellsize) * 5 * 5);
	for (k = 0; k < (int) (hight / cellsize) * atoi(argv[5]); k++) {
		for (i = 1; i < nrows * atoi(argv[4]); i++) {
			for (j = 1; j < ncols * atoi(argv[3]); j++) {
				if ((ind_multipl[i * ncols * atoi(argv[3]) + j] != -1) &&
					(ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j] != -1) &&
					(ind_multipl[i * ncols * atoi(argv[3]) + j - 1] != -1) &&
					(ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1] != -1)) {
						fprintf(f, "%d %d %d %d %d\n", 4,
							k * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1],
							k * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j],
							k * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1]);
						fprintf(f, "%d %d %d %d %d\n", 4,
							k * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j - 1],
							k * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j],
							k * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j],
							(k + 1) * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j]);
						fprintf(f, "%d %d %d %d %d\n", 4,
							k * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j]);
						fprintf(f, "%d %d %d %d %d\n", 4,
							k * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j],
							(k + 1) * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j]);
						fprintf(f, "%d %d %d %d %d\n", 4,
							k * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j],
							k * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[(i - 1) * ncols * atoi(argv[3]) + j - 1],
							(k + 1) * n_points_multipl + ind_multipl[i * ncols * atoi(argv[3]) + j]);
						a++;
				}
			}
		}
	}
	fprintf(f, "CELL_TYPES %d\n", n_cells * atoi(argv[3]) * atoi(argv[4]) * atoi(argv[5]) * (int) (hight / cellsize) * 5);
	for (i = 0; i < n_cells * atoi(argv[3]) * atoi(argv[4]) * atoi(argv[5]) * (int) (hight / cellsize) * 5; i++) {
		fprintf(f, "%d\n", 10);
	}
	return 0;
}

int free_massives(int argc, char **argv)
{
	extern float hight;
	extern int ncols;
	extern int nrows;
	extern float cellsize;
	extern float *mass;
	extern int *ind;
	extern int *bl_cond;
	extern int n_points, n_cells;
	extern int n_points_multipl;
	extern int *ind_multipl;
	extern int count_ind_multipl;
	extern float interpolation, interpolation_poli;
	extern float *mass1;

	free(mass);
	free(bl_cond);
	free(ind);
	free(ind_multipl);
	free(mass1);
	return 0;
}

int main(int argc, char **argv)
{
	extern float hight;
	extern int ncols;
	extern int nrows;
	extern float cellsize;
	extern float *mass;
	extern int *ind;
	extern int *bl_cond;
	extern int n_points, n_cells;
	extern int n_points_multipl;
	extern int *ind_multipl;
	extern int count_ind_multipl;
	extern float interpolation, interpolation_poli;
	extern float *mass1;

	if (argc != 7) goto error;

	if ((float) atoi(argv[3]) != atof(argv[3])) goto error;
	if ((float) atoi(argv[4]) != atof(argv[4])) goto error;
	if ((float) atoi(argv[5]) != atof(argv[5])) goto error;

	if (read_asc_and_declare_variables(argc, argv) == 1) goto error;
	if (do_interpolation(argc, argv) == 1) goto error;
	if (print_vtk(argc, argv) == 1) goto error;
	if (free_massives(argc, argv) == 1) goto error;

	return 0;
error:
	printf("Error\n");
	printf("asc2vtk.out [map.asc] [hight] [x koef] [y koef] [z koef]\n");
	return 1;
}
