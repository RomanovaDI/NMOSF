#pragma once

class init
{
private:
	char InputFileName[100];
	char MapName[100];
	char RegionFileName[100];
	int valDebugLevel;
public:
	init(char [100]);
	void ReadInputFile(char [100]);
	int DebugLevel();
};
