#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Mesh.h"
#include "MathModel.h"
#include "Init.h"

void init::readInputFile()
{
	FILE *inputFile = fopen(inputFileName, "r");
	if (inputFile == NULL)
		DEBUGP(0, FILE_OPEN_ERR, inputFileName);
	char str[300];
	char key[100];
	while (fgets(str, 300, inputFile)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "MapFileName")) {
				if (! sscanf(str, "%s %s", key, mapFileName))
					DEBUGP(0, FILE_DATA_ERR, "error MapFileName in file ", inputFileName);
			} else if (! strcmp(key, "RegionFileName")) {
				if (! sscanf(str, "%s %s", key, regionFileName))
					DEBUGP(0, FILE_DATA_ERR, "error RegionFileName in file ", inputFileName);
			} else if (! strcmp(key, "Dimension")) {
				if (! sscanf(str, "%s %d", key, &dimension))
					DEBUGP(0, FILE_DATA_ERR, "error Dimension in file ", inputFileName);
				if ((dimension < 1) || (dimension > 3))
					DEBUGP(0, FILE_DATA_ERR, "error Dimension in file ", inputFileName);
			} else if (! strcmp(key, "MeshCellSize")) {
				double cellsize;
				if (! sscanf(str, "%s %lf", key, &cellsize))
					DEBUGP(0, FILE_DATA_ERR, "error MeshCellSize in file ", inputFileName);
				Mesh.setCellSize(cellsize);
			} else if (! strcmp(key, "MathModelFileName")) {
				if (! sscanf(str, "%s %s", key, mathModelFileName))
					DEBUGP(0, FILE_DATA_ERR, "error RegionFileName in file ", inputFileName);
			} else
				DEBUGP(0, FILE_DATA_ERR, "unknown tag", key , "in file", inputFileName);	
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", inputFileName, "in line :", str);
	}
	fclose(inputFile);
}

void init::readMesh()
{
	Mesh.readASCII(mapFileName, regionFileName);
}

void init::printVTK()
{
	Mesh.printVTK();
}

void init::printInfo()
{
	cout << inputFileName << endl;
	cout << mapFileName << endl;
	cout << regionFileName << endl;
}

void init::setInputFileName(char name[100])
{
	strcpy(inputFileName, name);
}

void init::readMathModel()
{
	MathModel.readMathModel(mathModelFileName, dimension);
}
