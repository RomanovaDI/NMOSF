#pragma once

class init
{
private:
	char inputFileName[100];
	char mapFileName[100];
	char regionFileName[100];
public:
	mesh Mesh;
	init() {}
	void readInputFile();
	void readMesh();
	void printVTK();
	void printInfo();
	void setInputFileName(char [100]);
};
