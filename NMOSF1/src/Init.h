#pragma once

class init : public error
{
private:
	char InputFileName[100];
	char MapFileName[100];
	char RegionFileName[100];
public:
	init() : error() {}
	void ReadInputFile();
	void PrintInfo();
	void SetInputFileName(char [100]);
};
