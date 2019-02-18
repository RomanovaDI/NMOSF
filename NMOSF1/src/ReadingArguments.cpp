#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Mesh.h"
#include "MathModel.h"
#include "Init.h"
#include "ReadingArguments.h"

void ReadArgs(int argc, char *argv[], init *I)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-i")) {
			I->setInputFileName(argv[++i]);
		} else
			DEBUGP(0, INPUT_DATA_ERR, "unknown key", argv[i]);
	}
}
