#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Init.h"

void init::ReadInputFile()
{
	FILE *InputFile = fopen(InputFileName, "r");
	if (InputFile == NULL)
		DebugP(0, FILE_OPEN_ERR + ": " + InputFileName);
	char str[300];
	char key[100];
	while (fgets(str, 300, InputFile)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "RegionFileName")) {
				char regionFileName[100];
				if (! sscanf(str, "%s %s", key, regionFileName))
					cout << "error regionFileName" << endl;//DebugCout(0, FILE_DATA_ERR + ": error RegionFileName in file " + InputFileName);
				strcpy(RegionFileName, regionFileName);
			} else
				cout << "unknown tag" << endl;//DebugCout(0, FILE_DATA_ERR + ": unknown tag " + key + " in file " + InputFileName);	
		} else
			cout << "error key" << endl;
	}
	fclose(InputFile);
}

void init::PrintInfo()
{
	cout << InputFileName << endl;
	cout << MapName << endl;
	cout << RegionFileName << endl;
	cout << valDebugLevel << endl;
}

void init::SetMapName(char *name)
{
	strcpy(name, MapName);
}
