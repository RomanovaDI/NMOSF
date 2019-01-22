#include <iostream>
using namespace std;
#include "Error.h"

mesh::mesh(int nx, int ny, int nz, double cellSize)
{
	meshNx = nx;
	meshNy = ny;
	meshNz = nz;
	meshCellSize = cellSize;
	if ((meshCellInd = (int *) malloc(sizeof(int) * meshNX * meshNy * meshNz)) == NULL)
		DebugCout(0, MEM_ERR);
	return;
}

mesh::Nx()
{
	return meshNx;
}

mesh::Ny()
{
	return meshNy;
}

mesh::Nz()
{
	return meshNz;
}

void mesh::setMeshCellInd(int i, int j, int k, int ind)
{
	meshCellInd[k * meshNx * meshNy + j * meshNx + i] = ind;
	return;
}

int mesh::getMeshCellInd(int i, int j, int k)
{
	return meshCellInd[k * meshNx * meshNy + j * meshNx + i];
}

void mesh::setMeshNumActiveCells(int n)
{
	meshNumActiveCells = n;
	return;
}
