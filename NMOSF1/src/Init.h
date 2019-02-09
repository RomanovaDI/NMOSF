#pragma once

class init
{
private:
	char InputFileName[100];
	char MapFileName[100];
	char RegionFileName[100];
public:
	mesh Mesh;
	init() {}
	void ReadInputFile();
	void PrintInfo();
	void SetInputFileName(char [100]);
};
