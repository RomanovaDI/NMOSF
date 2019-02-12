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
	for (int i = 0; i < ncols * nrows; i++)
		if (fscanf(f, "%lf", &mass[i]) == EOF)
			DEBUGP(0, FILE_DATA_ERR, mapName);
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
	struct map Map;
	Map.mass = mass;
	Map.ncols = ncols;
	Map.nrows = nrows;
	Map.cellSize = cellsize;
	Map.nodataValue = nodata_value;
	if (meshCellSize != cellsize) {
		struct map MapOut;
		MapOut.ncols = floor((double) Map.ncols * Map.cellSize / meshCellSize);
		MapOut.nrows = floor((double) Map.nrows * Map.cellSize / meshCellSize);
		double *mass_interpolation;
		if ((mass_interpolation = (int *) malloc(MapOut.ncols * MapOut.nrows * sizeof(int))) == NULL)
			DEBUGP(0, MEM_ERR);
		MapOut.mass = mass_interpolation;
		MapOut.cellSize = meshCellSize;
		MapOut.nodataValue = nodata_value;
		interpolate(&Map, &MapOut);
		Map = MapOut;
		free(mass);
	}
	double max = 0, min = 0;
	for (int i = 0; i < Map.nrows; i++)
		for (int j = 0; j < Map.ncols; j++) {
			if ((max == 0) && (min == 0) && (Map.mass[i * Map.ncols + j] != Map.nodataValue))
				max = min = Map.mass[i * Map.ncols + j];
			max = (max < Map.mass[i * Map.ncols + j]) ? Map.mass[i * Map.ncols + j] : max;
			min = (min > Map.mass[i * Map.ncols + j]) ? Map.mass[i * Map.ncols + j] : min;
		}
	meshNx = Map.nrows;
	meshNy = Map.ncols;
	double hight = 20;
	meshNz = ceil ((max - min + hight) / Map.cellSize);
	if (meshCellInd = (int *) malloc(meshNx * meshNy * meshNz * sizeof(int)) == NULL)
		DebugP(0, MEM_ERR);
	for (int k = 0, meshNumActiveCells = 0; k < meshNz; k++)
		for (int i = 0; i < meshNx; i++)
			for (int j = 0; j < meshNy; j++)
				meshCellInd[k * meshNx * meshNy + i * meshNy + j] =
					((Map.mass[i * Map.ncols + j] != Map.nodataValue) && (k * meshCellSize >= Map.mass[i * Map.ncols + j]) && (k * meshCellSize < Map.mass[i * Map.ncols + j] + hight)) ?
					meshNumActiveCells++ : -1;
}

void mesh::interpolate(struct map *MapIn, struct map *MapOut)
{
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
				double a = MapIn.mass[i * MapIn.ncols + j];
				double d = (c[j] - c[j - 1]) / MapIn.cellsize;
				double b = (MapIn.mass[i * MapIn.ncols + j] - MapIn.mass[i * MapIn.ncols + j - 1]) / MapIn.cellsize + MapIn.cellsize * (2 * c[j] + c[j - 1]) / 6;
				int ind = ceil((double) (j - 1) * MapIn.cellSize / MapOut.cellSize);
				while (ind * MapOut.cellSize <= j * MapIn.cellSize) {
					mass[i * ncols + ind] = a + b * (ind * MapOut.cellSize - j * MapIn.cellSize) + c[j] *
						pow(ind * MapOut.cellSize - j * MapIn.cellSize, 2) + d * pow(ind * MapOut.cellSize - j * MapIn.cellSize, 3);
					ind++;
				}
			}
		}
	nrows = MapOut.nrows;
	double *mass1;
	if ((mass1 = (double *) malloc(ncols * nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	for (i = 0; i < ncols * nrows; i++)
		mass1[i] = MapOut.nodataValue;
	if ((c = (double *) realloc(c, MapIn.nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if ((e = (double *) realloc(e, MapIn.nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if ((f = (double *) realloc(f, MapIn.nrows * sizeof(double))) == NULL)
		DEBUGP(0, MEM_ERR);
	if (MapIn.nrows > 2)
		for (int j = 0; j < MapOut.ncols; j++) {
			flag = 0;
			ind_start = 0;
			ind_finish = MapIn.nrows;
			for (int i = 0; i < MapIn.nrows; i++) {
				if ((mass[i * MapOut.ncols + j] != MapIn.nodataValue) && (flag++ == 0))
					ind_start = i;
				if ((mass[i * MapOut.ncols + j] == MapIn.nodataValue) && (flag++ == 1))
					ind_finish = i;
			}
			memset(e, 0, MapIn.ncols * sizeof(double));
			memset(f, 0, MapIn.ncols * sizeof(double));
			for (int i = ind_start + 3; i < ind_finish - 1; i++) {
				e[i] = - 1. / (4 + e[i - 1]);
				f[i] = (3 * (mass[(i + 1) * MapOut.ncols + j] - 2 * mass[i * MapOut.ncols + j] + mass[(i - 1) * MapOut.ncols + j]) / (MapIn.cellSize * MapIn.cellSize) - f[i - 1]) / (4 + e[i - 1]);
			}
			memset(c, 0, MapIn.nrows * sizeof(double));
			i = ind_finish - 1;
			c[i] = (3 * (mass[(i + 1) * MapOut.ncols + j] - 2 * mass[i * MapOut.ncols + j] + mass[(i - 1) * MapOut.ncols + j]) / (MapIn.cellSize * MapIn.cellSize) - f[i]) / (4 + e[i]);
			for (int i = ind_finish - 2; i > ind_start + 1; i--)
				c[i] = e[i + 1] * c[i + 1] + f[i + 1];
			for (int i = ind_start + 1; i < ind_finish; i++) {
				double a = mass[i * MapOut.ncols + j];
				double d = (c[i] - c[i - 1]) / MapIn.cellsize;
				double b = (mass[i * MapOut.ncols + j] - mass[(i - 1) * MapOut.ncols + j]) / MapIn.cellsize + MapIn.cellsize * (2 * c[i] + c[i - 1]) / 6;
				int ind = ceil((double) (i - 1) * MapIn.cellSize / MapOut.cellSize);
				while (ind * MapOut.cellSize <= i * MapIn.cellSize) {
					mass1[ind * ncols + i] = a + b * (ind * MapOut.cellSize - i * MapIn.cellSize) + c[i] * pow(ind * MapOut.cellSize - i * MapIn.cellSize, 2) +
						d * pow(ind * MapOut.cellSize - i * MapIn.cellSize, 3);
					ind++;
				}
			}
		}
	for (int i = 0; i < MapOut.nrows; i++)
		for (int j = 0; j < MapOut.ncols; j++)
			MapOut.mass[i * MapOut.ncols + j] = mass1[i * MapOut.ncols + j];
	free(c);
	free(e);
	free(f);
	free(mass);
	free(mass1);
}

void mesh::printVTK()
{
	FILE *f = fopen("map.vtk", "w");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, "map.vtk");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Map of slope, -1 is unused cells\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d %d %d\n", meshNx, meshNy, meshNz);
}

void mesh::setCellSize(double cellsize)
{
	meshCellSize = cellsize;
}

double mesh::getCellSize()
{
	return meshCellSize;
}
