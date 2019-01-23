#include <iostream>
#include <stdio.h>
using namespace std;
//#include "MathModel.h"
#include "Error.h"
#include "Mesh.h"

int main(int argc, char* argv[])
{
	if (argc == 1)
		DebugCout(0, INIT_DATA_ERR);
	char initData[100] = "relief_22.asc";
	mesh myMesh(initData);
    return 0;
}
