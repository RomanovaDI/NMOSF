#include <iostream>
#include <stdio.h>
using namespace std;
//#include "MathModel.h"
#include "Error.h"
#include "Init.h"
//#include "Mesh.h"

int main(int argc, char* argv[])
{
	if (argc != 2)
		DebugCout(0, INIT_DATA_ERR);
	init Init(argv[1]);
	//mesh myMesh(initData);
    return 0;
}
