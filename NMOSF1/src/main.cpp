#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Init.h"
#include "ReadingArguments.h"

int main(int argc, char* argv[])
{
	cout << "starting" << endl;
	init II;
	init *I = &II;
	ReadArgs(argc, argv, I);
	I->ReadInputFile();
	cout << "input file was read" << endl;
	int i = 4;
	I->DebugP(0, MEM_ERR, "Debugging", i);
	I->PrintInfo();
	//mesh myMesh(initData);
    return 0;
}
