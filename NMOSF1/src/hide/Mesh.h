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
	int meshTmpInd;
public:
	mesh(char[100]);
	int Nx();
	int Ny();
	int Nz();
	void setMeshCellInd(int, int, int);
	int getMeshCellInd(int, int, int);
	void setMeshNumActiveCells(int);
	int readASCII(char[100]);
};