#pragma once

struct map
{
	double *mass;
	int ncols, nrows;
	double cellSize;
	double nodataValue;
}

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
	double cellSize;
	void interpolate(struct map *, struct map *);
public:
	mesh(double);
	int Nx();
	int Ny();
	int Nz();
	//void setMeshCellInd(int, int, int);
	//int getMeshCellInd(int, int, int);
	//void setMeshNumActiveCells(int);
	void readASCII(char[100], char[100]);
	void setCellSize(double);
	double getCellSize();
};
