#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Mesh.h"
#include "Init.h"

void init::ReadInputFile()
{
	FILE *InputFile = fopen(InputFileName, "r");
	if (InputFile == NULL)
		DEBUGP(0, FILE_OPEN_ERR, InputFileName);
	char str[300];
	char key[100];
	while (fgets(str, 300, InputFile)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "MapFileName")) {
				char mapFileName[100];
				if (! sscanf(str, "%s %s", key, mapFileName))
					DEBUGP(0, FILE_DATA_ERR, "error MapFileName in file ", InputFileName);
				strcpy(MapFileName, mapFileName);
			} else if (! strcmp(key, "RegionFileName")) {
				char regionFileName[100];
				if (! sscanf(str, "%s %s", key, regionFileName))
					DEBUGP(0, FILE_DATA_ERR, "error RegionFileName in file ", InputFileName);
				strcpy(RegionFileName, regionFileName);
			} else if (! strcmp(key, "MeshCellSize")) {
				double cellsize;
				if (! sscanf(str, "%s %lf", key, &cellsize))
					DEBUGP(0, FILE_DATA_ERR, "error MeshCellSize in file ", InputFileName);
				Mesh.setCellSize(cellsize);
			} else
				DEBUGP(0, FILE_DATA_ERR, "unknown tag", key , "in file", InputFileName);	
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", InputFileName, "in line :", str);
	}
	fclose(InputFile);
}

void init::PrintInfo()
{
	cout << InputFileName << endl;
	cout << MapFileName << endl;
	cout << RegionFileName << endl;
}

void init::SetInputFileName(char name[100])
{
	strcpy(InputFileName, name);
}
