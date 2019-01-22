#pragma once

class mesh
{
private:
	int meshNx;
	int meshNy;
	int meshNz;
	int *meshCellInd;
	double meshCellSize;
	int meshNumActiveCells;
public:
	mesh(int, int, int, double);
	int Nx();
	int Ny();
	int Nz();
	void setMeshCellInd(int, int, int, int);
	int getMeshCellInd(int, int, int);
	void setMeshNumActiveCells(int);
}
