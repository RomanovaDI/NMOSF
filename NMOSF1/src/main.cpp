#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Init.h"

int main(int argc, char* argv[])
{
	cout << "starting" << endl;
	init II;
	init *I = &II;
	cout << "Correct numbeg of arguments" << endl;
	cout << "initializing init" << endl;
	I->ReadInputFile();
	cout << "input file was read" << endl;
	I->PrintInfo();
	//mesh myMesh(initData);
    return 0;
}
