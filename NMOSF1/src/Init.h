#pragma once

class error;

class init : public error
{
private:
	char InputFileName[100];
	char MapName[100];
	char RegionFileName[100];
public:
	init() : error() {}
	void ReadInputFile();
	void PrintInfo();
	void SetMapName(char *);
};
