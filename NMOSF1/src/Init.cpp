#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Init.h"

init::init(char name[100])
{
	ReadInputFile(name);
}

void init::ReadInputFile(char name[100])
{
	strcpy(InputFileName, name);
	FILE *InputFile = fopen(InputFileName, "r");
	if (InputFile == NULL)
		DebugCout(0, FILE_OPEN_ERR << ": " << InputFileName);
	char str[300];
	char key[100];
	while (fgets(str, 300, InputFile)) {
		if (str[0] == "#")
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "DEBUG")) {
				int debugLevel;
				if (! sscanf(str, "%d", &debugLevel))
					DebugCout(0, FILE_DATA_ERR << ": error debug level in file " << InputFileName);
				valDebugLevel = debugLevel;
			} else if (! strcmp(key, "MapName")) {
				char mapName[100];
				if (! sscanf(str, "%s", mapName))
					DebugCout(0, FILE_DATA_ERR << ": error MapName in file " << InputFileName);
				strcpy(MapName, mapName);
			} else if (! strcmp(key, "RegionFileName")) {
				char regionFileName[100];
				if (! sscanf(str, "%s", regionFileName))
					DebugCout(0, FILE_DATA_ERR << ": error RegionFileName in file " << InputFileName);
				strcpy(RegionFileName, regionFileName);
			} else
				DebugCout(0, FILE_DATA_ERR << ": unknown tag " << key << " in file " << InputFileName);	
		}
	}
	fclose(InputFile);
}

int init::DeugLevel()
{
	return valDebugLevel;
}
