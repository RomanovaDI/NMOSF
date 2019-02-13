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
#include "ReadingArguments.h"

int main(int argc, char* argv[])
{
	cout << "starting" << endl;
	init II;
	init *I = &II;
	cout << "variables were created" << endl;
	ReadArgs(argc, argv, I);
	cout << "args were red" << endl;
	I->readInputFile();
	cout << "input file was read" << endl;
	int i = 4;
	I->printInfo();
	I->readMesh();
	I->printVTK();
	//mesh myMesh(initData);
    return 0;
}
