#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>
using namespace std;
#include "Error.h"
#include "Mesh.h"

mesh::mesh(double cellsize)
{
	setCellSize(cellsize);
}

int mesh::Nx()
{
	return meshNx;
}

int mesh::Ny()
{
	return meshNy;
}

int mesh::Nz()
{
	return meshNz;
}

/*
void mesh::setMeshCellInd(int i, int j, int k)
{
	meshCellInd[k * meshNx * meshNy + j * meshNx + i] = meshTmpInd++;
}

int mesh::getMeshCellInd(int i, int j, int k)
{
	return meshCellInd[k * meshNx * meshNy + j * meshNx + i];
}

void mesh::setMeshNumActiveCells(int n)
{
	meshNumActiveCells = n;
}
*/

void mesh::readASCII(char mapName[100], char regionName[100])
{
	FILE *f = fopen(mapName,"r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, mapName);
	FILE *f1 = fopen("map.txt", "w");
	while ((int i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);
	f = fopen(regionName,"r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, regionName);
	f1 = fopen("regions_map.txt", "w");
	while ((int i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);
	f = fopen("map.txt", "r");
	f1 = fopen("regions_map.txt", "r");
	char str[20];
	int ncols;
	if (fscanf(f, "%s %d", str, &ncols) == EOF)
		DEBUGP(0, FILE_DATA_ERR, mapName);
	if (strcmp(str, "ncols") != 0) {
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"ncols\" tag");
	int nrows;
	if (fscanf(f, "%s %d", str, &nrows) == EOF)
		DEBUGP(0, FILE_DATA_ERR, mapName);
	if (strcmp(str, "nrows") != 0) {
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"nrows\" tag");
	double xllcorner;
	if (fscanf(f, "%s %lf", str, &xllcorner) == EOF)
		DEBUGP(0, FILE_DATA_ERR, mapName);
	if (strcmp(str, "xllcorner") != 0) {
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"xllcorner\" tag");
	double yllcorner;
	if (fscanf(f, "%s %lf", str, &yllcorner) == EOF)
		DEBUGP(0, FILE_DATA_ERR, mapName);
	if (strcmp(str, "yllcorner") != 0) {
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"yllcorner\" tag");
	double cellsize;
	if (fscanf(f, "%s %lf", str, &cellsize) == EOF)
		DEBUGP(0, FILE_DATA_ERR, mapName);
	if (strcmp(str, "cellsize") != 0)
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"cellsize\" tag");
	int nodata_value;
	if (fscanf(f, "%s %lf", str, &nodata_value) == EOF)
		DEBUGP(0, FILE_DATA_ERR << " " << mapName);
	if (strcmp(str, "NODATA_value") != 0)
		DEBUGP(0, FILE_DATA_ERR, mapName, ": no \"NODATA_value\" tag");
	int ncols1;
	if (fscanf(f1, "%s %d", str, &ncols1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "ncols") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"ncols\" tag");
	int nrows1;
	if (fscanf(f1, "%s %d", str, &nrows1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "nrows") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"nrows\" tag");
	double xllcorner1;
	if (fscanf(f1, "%s %lf", str, &xllcorner1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "xllcorner") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"xllcorner\" tag");
	double yllcorner1;
	if (fscanf(f1, "%s %lf", str, &yllcorner1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "yllcorner") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"yllcorner\" tag");
	double cellsize1;
	if (fscanf(f1, "%s %lf", str, &cellsize1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "cellsize") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"cellsize\" tag");
	double nodata_value1;
	if (fscanf(f1, "%s %lf", str, &nodata_value1) == EOF)
		DEBUGP(0, FILE_DATA_ERR, regionName);
	if (strcmp(str, "NODATA_value") != 0)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": no \"NODATA_value\" tag");
	if (cellsize != cellsize1)
		DEBUGP(0, FILE_DATA_ERR, regionName, ": Cellsize in both maps need to be the same");
	double *mass;
	if ((mass = (double *) malloc(ncols * nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	memset((void *) mass, 0, ncols * nrows * sizeof(double));
	/*
	if ((I->ind = (int *) malloc(I->ncols * I->nrows * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->ind, 0, I->ncols * I->nrows * sizeof(int));
	if ((I->bl_cond = (int *) malloc((I->ncols - 1) * (I->nrows - 1) * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
	memset((void *) I->bl_cond, 0, (I->ncols - 1) * (I->nrows - 1) * sizeof(int));
	*/
	for (int i = 0; i < ncols * nrows; i++) {
		if (fscanf(f, "%lf", &mass[i]) == EOF)
			DEBUGP(0, FILE_DATA_ERR, mapName);
		/*
		I->ind[i] = -1;
		if (i < (I->ncols - 1) * (I->nrows - 1))
			I->bl_cond[i] = -1;
		*/
	}
	fclose(f);
	remove("map.txt");
	double *mass_tmp;
	if ((mass_tmp = (double *) malloc(ncols1 * nrows1 * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	memset((void *) mass_tmp, 0, ncols1 * nrows1 * sizeof(double));
	for (int i = 0; i < ncols1 * nrows1; i++)
		if (fscanf(f1, "%lf", &mass_tmp[i]) == EOF)
			DEBUGP(0, FILE_DATA_ERR, regionName);
	remove("regions_map.txt");
	if ((snow_region = (int *) malloc(ncols * nrows * sizeof(int))) == NULL)
		DEBUGP(0, MEM_ERR);
	memset((void *) snow_region, 0, ncols * nrows * sizeof(int));
	double xlucorner = xllcorner - cellsize * nrows;
	double ylucorner = yllcorner;
	double xlucorner1 = xllcorner1 - cellsize1 * nrows1;
	double ylucorner1 = yllcorner1;
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			if ((xlucorner1 - xlucorner >= 0) && (ylucorner -ylucorner1 >= 0) &&
				(j * cellsize + ylucorner >= ylucorner1) &&
				(i * cellsize + xlucorner >= xlucorner1) &&
				(j * cellsize + ylucorner <= ylucorner1 + ncols1 * cellsize) &&
				(i * cellsize + xlucorner <= xlucorner1 + nrows1 * cellsize)) {
					int tmp_ind = (i - (int) ((xlucorner1 - xlucorner) / cellsize)) * ncols1 + j - (int) ((ylucorner1 - ylucorner) / cellsize);
					if (mass_tmp[tmp_ind] == 1)
						snow_region[i * ncols + j] = 1;
					else if (mass_tmp[tmp_ind] == 0)
						snow_region[i * ncols + j] = 0;
					else
						snow_region[i * ncols + j] = -1;
			} else {
					snow_region[i * ncols + j] = -1;
			}
		}
	}
	free(mass_tmp);
	for (int i = 0; i < nrows; i++) {
		int flag = 0;
		for (int j = 0; j < ncols; j++) {
			if ((flag % 2 == 0) && (mass[i * ncols + j] != nodata_value))
				flag++;
			if ((flag % 2 == 1) && (mass[i * ncols + j] == nodata_value))
				flag++;
		}
		if (flag > 2)
			DEBUGP(0, MAP_FILE_ERR, "holes in map in row", i);
	}
	for (int j = 0; j < ncols; j++) {
		int flag = 0;
		for (int i = 0; i < nrows; i++) {
			if ((flag % 2 == 0) && (mass[i * ncols + j] != nodata_value))
				flag++;
			if ((flag % 2 == 1) && (mass[i * ncols + j] == nodata_value))
				flag++;
		}
		if (flag > 2)
			DEBUGP(0, MAP_FILE_ERR, "holes in map in column", j);
	}
	struct map MapIn;
	MapIn.mass = mass;
	MapIn.ncols = ncols;
	MapIn.nrows = nrows;
	MapIn.cellSize = cellsize;
	MapIn.nodataValue = nodata_value;
	struct map MapOut;
	MapOut.ncols = floor((double) Map.ncols * Map.cellSize / cellSize);
	MapOut.nrows = floor((double) Map.nrows * Map.cellSize / cellSize);
	double *mass_interpolation;
	if ((mass_interpolation = (int *) malloc(MapOut.ncols * MapOut.nrows * sizeof(int))) == NULL)
		DEBUGP(0, MEM_ERR);
	MapOut.mass = mass_interpolation;
	MapOut.cellSize = cellSize;
	MapOut.nodataValue = nodata_value;
	if (cellSize != cellsize)
		interpolate(&MapIn, &MapOut);
	else
		MapOut = MapIn;
}

void mesh::interpolate(struct map *MapIn, struct map *MapOut)
{
	double a, b, d;
	int ncols = MapOut.ncols;
	int nrows = MapIn.nrows;
	double *mass;
	if ((mass = (double *) malloc(ncols * nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	for (i = 0; i < ncols * nrows; i++)
		mass[i] = MapOut.nodataValue;
	double *c;
	if ((c = (double *) malloc(MapIn.ncols * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	double *e;
	if ((e = (double *) malloc(MapIn.ncols * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	double *f;
	if ((f = (double *) malloc(MapIn.ncols * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	int ind_start, ind_finish, flag;
	if (MapIn.ncols > 2)
		for (int i = 0; i < MapIn.nrows; i++) {
			flag = 0;
			ind_start = 0;
			ind_finish = MapIn.ncols;
			for (int j = 0; j < MapIn.ncols; j++) {
				if ((MapIn.mass[i * MapIn.ncols + j] != MapIn.nodataValue) && (flag++ == 0))
					ind_start = j;
				if ((MapIn.mass[i * MapIn.ncols + j] == MapIn.nodataValue) && (flag++ == 1))
					ind_finish = j;
			}
			memset(e, 0, MapIn.ncols * sizeof(double));
			memset(f, 0, MapIn.ncols * sizeof(double));
			for (int j = ind_start + 3; j < ind_finish - 1; j++) {
				e[j] = - 1. / (4 + e[j - 1]);
				f[j] = (3 * (MapIn.mass[i * MapIn.ncols + j + 1] - 2 * MapIn.mass[i * MapIn.ncols + j] + MapIn.mass[i * MapIn.ncols + j - 1]) / (MapIn.cellSize * MapIn.cellSize) - f[j - 1]) / (4 + e[j - 1]);
			}
			memset(c, 0, MapIn.ncols * sizeof(double));
			j = ind_finish - 1;
			c[j] = (3 * (MapIn.mass[i * MapIn.ncols + j + 1] - 2 * MapIn.mass[i * MapIn.ncols + j] + MapIn.mass[i * MapIn.ncols + j - 1]) / (MapIn.cellSize * Map.IncellSize) - f[j]) / (4 + e[j]);
			for (int j = ind_finish - 2; j > ind_start + 1; j--)
				c[j] = e[j + 1] * c[j + 1] + f[j + 1];
			for (int j = ind_start + 1; j < ind_finish; j++) {
				a = MapIn.mass[i * MapIn.ncols + j];
				d = (c[j] - c[j - 1]) / MapIn.cellsize;
				b = (MapIn.mass[i * MapIn.ncols + j] - MapIn.mass[i * MapIn.ncols + j - 1]) / MapIn.cellsize + MapIn.cellsize * (2 * c[j] + c[j - 1]) / 6;
				int ind = ceil((double) (j - 1) * MapIn.cellSize / MapOut.cellSize);
				while (ind * MapOut.cellSize <= j * MapIn.cellSize) {
					mass[i * ncols + ind] = a + b * (ind * MapOut.cellSize - j * MapIn.cellSize) + c[j] * pow(ind * MapOut.cellSize - j * MapIn.cellSize, 2) + d * pow(ind * MapOut.cellSize - j * MapIn.cellSize, 3);
					ind++;
				}
			}
		}
	///////////////////////////////////////////////////////////////////////////////////
	Map.mass = mass;
	nrows = floor((double) Map.nrows * Map.cellSize / cellSize);
	if ((c = (double *) realloc(c, I->nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if ((e = (double *) realloc(e, I->nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if ((f = (double *) realloc(f, I->nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if (Map.nrows > 2)
		for (int j = 0; j < Map.ncols; j++) {
			flag = 0;
			ind_start = 0;
			ind_finish = Map.nrows;
			for (int i = 0; i < Map.nrows; i++) {
				if ((Map.mass[i * Map.ncols + j] != Map.nodataValue) && (flag++ == 0))
					ind_start = i;
				if ((Map.mass[i * Map.ncols + j] == Map.nodataValue) && (flag++ == 1))
					ind_finish = i;
			}
			memset(e, 0, Map.ncols * sizeof(double));
			memset(f, 0, Map.ncols * sizeof(double));
			for (int i = ind_start + 3; i < ind_finish - 1; i++) {
				e[i] = - 1. / (4 + e[i - 1]);
				f[i] = (3 * (Map.mass[(i + 1) * Map.ncols + j] - 2 * Map.mass[i * Map.ncols + j] + Map.mass[(i - 1) * Map.ncols + j]) / (Map.cellSize * Map.cellSize) - f[i - 1]) / (4 + e[i - 1]);
			}
			memset(c, 0, Map.nrows * sizeof(double));
			i = ind_finish - 1;
			c[i] = (3 * (Map.mass[(i + 1) * Map.ncols + j] - 2 * Map.mass[i * Map.ncols + j] + Map.mass[(i - 1) * Map.ncols + j]) / (Map.cellSize * Map.cellSize) - f[i]) / (4 + e[i]);
			for (int i = ind_finish - 2; i > ind_start + 1; i--)
				c[i] = e[i + 1] * c[i + 1] + f[i + 1];
			for (int i = ind_start + 1; i < ind_finish; i++) {
				a = Map.mass[i * I->ncols + j];
				d = (c[i] - c[i - 1]) / Map.cellsize;
				b = (Map.mass[i * Map.ncols + j] - Map.mass[(i - 1) * Map.ncols + j]) / Map.cellsize + Map.cellsize * (2 * c[i] + c[i - 1]) / 6;
				int ind = ceil((double) (i - 1) * Map.cellSize / cellSize);
				while (ind * cellSize <= i * Map.cellSize) {
					mass[ind * ncols + i] = a + b * (ind * cellSize - i * Map.cellSize) + c[i] * pow(ind * cellSize - i * Map.cellSize, 2) + d * pow(ind * cellSize - i * Map.cellSize, 3);
					ind++;
				}
			}
		}

	if ((I->nrows > 2) && (I->kx > 1)) {
		for (j = 0; j < I->ncols; j++) {
			for (l = 0; l < I->ky; l++) {
				flag = 0;
				ind_start = 0;
				ind_finish = I->nrows;
				for (i = 0; i < I->nrows; i++) {
					if (((I->ind[i * I->ncols + j] != -1) && (flag == 0) && (l == 0)) ||
						(((I->ind[i * I->ncols + j] != -1) && ((j + 1 == I->ncols) || (I->ind[i * I->ncols + j + 1] != -1))) &&
						 	(flag == 0) && (l != 0))) {
								ind_start = i;
								flag = 1;
					}
					if (((I->ind[i * I->ncols + j] == -1) && (flag == 1) && (l == 0)) ||
						(((I->ind[i * I->ncols + j] == -1) || ((j + 1 == I->ncols) || (I->ind[i * I->ncols + j + 1] == -1))) &&
						 	(flag == 1) && (l != 0))) {
								ind_finish = i;
								flag = 2;
					}
				}
				e[ind_start + 2] = f[ind_start + 2] = 0;
				for (i = ind_start + 3; i < ind_finish - 1; i++) {
					e[i] = - 1. / (4 * e[i - 1] + 1);
					f[i] = (6 * (I->mass1[(i + 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
						2 * I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] +
						I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
						(I->cellsize * I->cellsize) - 4 * f[i - 1]) / (4 * e[i - 1] + 1);
				}
				i = ind_finish - 2;
				c[i] = (6 * (I->mass1[(i + 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
					2 * I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] +
					I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
					(I->cellsize * I->cellsize) - 4 * f[i]) / (1 + 4 * e[i]);
				for (i = ind_finish - 3; i > ind_start + 1; i--) {
					c[i] = e[i + 1] * c[i + 1] + f[i + 1];
				}
				c[ind_start + 1] = c[ind_finish - 1] = c[ind_start] = 0;
				for (i = ind_start + 1; i < ind_finish; i++) {
					a = I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l];
					d = (c[i] - c[i - 1]) / I->cellsize;
					b = (I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] -
						I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l]) /
						I->cellsize + I->cellsize * (2 * c[i] + c[i - 1]) / 6;
					//printf("j = %d, l = %d, i = %d, a = %f, b = %f, c = %f, d = %f\n", j, l, i, a, b, c[i], d);
					if ((I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] != I->nodata_value) &&
						(I->mass1[i * ((I->ncols - 1) * I->ky + 1) * I->kx + j * I->ky + l] != I->nodata_value)) {
							for (k = 0; k < I->kx; k++) {
								I->mass1[(i - 1) * ((I->ncols - 1) * I->ky + 1) * I->kx + k * ((I->ncols - 1) * I->ky + 1) + j * I->ky + l] =
									a + b * (k * I->cellsize / (double) I->kx - I->cellsize) +
									c[i] * pow(k * I->cellsize / (double) I->kx - I->cellsize, 2) / 2 +
									d * pow(k * I->cellsize / (double) I->kx - I->cellsize, 3) / 6;
							}
					}
				}
				memset((void *) c, 0, I->nrows * sizeof(double));
				memset((void *) e, 0, I->nrows * sizeof(double));
				memset((void *) f, 0, I->nrows * sizeof(double));
				if (j == I->ncols - 1) break;
			}
		}
#if DEBUG
		printf("Interpolation along the x axes have been done in process %d\n", I->my_rank);
#endif
	}

	free(c);
	free(e);
	free(f);

	count_ind_multipl = 0;
	for (i = 0; i < (I->nrows - 1) * I->kx + 1; i++) {
		for (j = 0; j < (I->ncols - 1) * I->ky + 1; j++) {
			if (I->mass1[i * ((I->ncols - 1) * I->ky + 1) + j] != I->nodata_value) {
				I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] = count_ind_multipl;
				count_ind_multipl++;
			} else {
				I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] = -1;
			}
		}
	}

	if ((I->ind_cell_multipl = (int *) malloc((I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky * sizeof(int))) == NULL) {
		printf("Memory error in function %s in process %d\n", __func__, I->my_rank);
		return 1;
	}
//	memset((void *) I->ind_cell_multipl, 0, (I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky * sizeof(int));

	for (i = 0; i < (I->nrows - 1) * I->kx * (I->ncols - 1) * I->ky; i++) {
		I->ind_cell_multipl[i] = -1;
	}

	I->n_cells_multipl = 0;
	for (i = 1; i < (I->nrows - 1) * I->kx + 1; i++) {
		for (j = 1; j < (I->ncols - 1) * I->ky + 1; j++) {
			if ((I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
				(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j] != -1) &&
				(I->ind_multipl[i * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1) &&
				(I->ind_multipl[(i - 1) * ((I->ncols - 1) * I->ky + 1) + j - 1] != -1)) {
					I->ind_cell_multipl[(i - 1) * (I->ncols - 1) * I->ky + j - 1] = I->n_cells_multipl;
					(I->n_cells_multipl)++;
			}
		}
	}

	I->dx[0] /= (double) I->kx;
	I->dx[1] /= (double) I->ky;
	I->dx[2] /= (double) I->kz;
}

void mesh::setCellSize(double cellsize)
{
	cellSize = cellsize;
}

double mesh::getCellSize()
{
	return cellSize;
}
