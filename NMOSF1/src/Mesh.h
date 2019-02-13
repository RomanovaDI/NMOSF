#pragma once

struct map
{
	double *mass;
	int ncols, nrows;
	double cellSize;
	double nodataValue;
};

class mesh
{
private:
	int meshNx;
	int meshNy;
	int meshNz;
	int *meshCellInd;
	double meshCellSize;
	int meshNumActiveCells;
	void interpolate(struct map *, struct map *);
public:
	mesh() {}
	int Nx();
	int Ny();
	int Nz();
	void readASCII(char[100], char[100]);
	void setCellSize(double);
	double getCellSize();
	void printVTK();
};
