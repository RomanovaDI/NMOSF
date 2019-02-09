#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Mesh.h"

//mesh::mesh(char fileName[100])
//{
//	mesh::readASCII(fileName);
/*	meshNx = nx;
	meshNy = ny;
	meshNz = nz;
	meshCellSize = cellSize;
	if ((meshCellInd = (int *) malloc(sizeof(int) * meshNx * meshNy * meshNz)) == NULL)
		DebugCout(0, MEM_ERR);
	memset(meshCellInd, 0, sizeof(int) * meshNx * meshNy * meshNz);
	meshTmpInd = 0;*/
//}

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
	int i, j;
	while ((i = getc(f)) != EOF) {
		if (i == ',') i = '.';
		putc(i, f1);
	}
	fclose(f1);
	fclose(f);
	f = fopen(regionName,"r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, regionName);
	f1 = fopen("regions_map.txt", "w");
	while ((i = getc(f)) != EOF) {
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
	for (i = 0; i < ncols * nrows; i++) {
		if (err = fscanf(f, "%lf", &mass[i]) == EOF)
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
	for (i = 0; i < ncols1 * nrows1; i++)
		if (err = fscanf(f1, "%lf", &mass_tmp[i]) == EOF)
			DEBUGP(0, FILE_DATA_ERR, regionName);
	remove("regions_map.txt");
	if ((snow_region = (int *) malloc(ncols * nrows * sizeof(int))) == NULL)
		DEBUGP(0, MEM_ERR);
	memset((void *) snow_region, 0, ncols * nrows * sizeof(int));
	double xlucorner = xllcorner - cellsize * nrows;
	double ylucorner = yllcorner;
	double xlucorner1 = xllcorner1 - cellsize1 * nrows1;
	double ylucorner1 = yllcorner1;
	for (i = 0; i < I->nrows; i++) {
		for (j = 0; j < I->ncols; j++) {
			if ((xlucorner1 - xlucorner >= 0) && (ylucorner -ylucorner1 >= 0) &&
				(j * I->cellsize + ylucorner >= ylucorner1) &&
				(i * I->cellsize + xlucorner >= xlucorner1) &&
				(j * I->cellsize + ylucorner <= ylucorner1 + ncols1 * I->cellsize) &&
				(i * I->cellsize + xlucorner <= xlucorner1 + nrows1 * I->cellsize)) {
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
	I->n_points = 0;
	I->n_cells = 0;
	for (i = 0; i < I->nrows; i++) {
		for (j = 0; j < I->ncols; j++) {
			if ((i - 1 >= 0) && (j + 1 < I->ncols) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j + 1] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j + 1] != I->nodata_value)) {
					I->ind[i * I->ncols +j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((j + 1 < I->ncols) && (i + 1 < I->nrows) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j + 1] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j + 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((i + 1 < I->nrows) && (j - 1 >= 0) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j - 1] != I->nodata_value) &&
				(I->mass[(i + 1) * I->ncols + j - 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(I->mass[i * I->ncols + j] != I->nodata_value) &&
				(I->mass[i * I->ncols + j - 1] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j] != I->nodata_value) &&
				(I->mass[(i - 1) * I->ncols + j - 1] != I->nodata_value)) {
					I->ind[i * I->ncols + j] = I->n_points;
					I->n_points++;
					goto end_points;
			}
			end_points:
			if ((j - 1 >= 0) && (i - 1 >= 0) &&
				(I->ind[i * I->ncols + j] != -1) &&
				(I->ind[i * I->ncols + j - 1] != -1) &&
				(I->ind[(i - 1) * I->ncols + j - 1] != -1) &&
				(I->ind[(i - 1) * I->ncols + j] != -1)) {
					I->bl_cond[(i - 1) * (I->ncols - 1) + j - 1] = I->n_cells;
					I->n_cells++;
			}
			if (I->ind[i * I->ncols + j] == -1)
				I->mass[i * I->ncols + j] = I->nodata_value;
		}
	}
	I->gl_nx = I->nx = (I->nrows - 1) * I->kx;
	I->gl_ny = I->ny = (I->ncols - 1) * I->ky;
	I->gl_nz = I->nz = (int) (I->hight / I->cellsize) * I->kz;
	I->dx[0] = I->dx[1] = I->dx[2] = I->cellsize;
}


void mesh::setCellSize(double cellsize)
{
	cellSize = cellsize;
}

double mesh::getCellSize()
{
	return cellSize;
}
