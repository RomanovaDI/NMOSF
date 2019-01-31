#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
using namespace std;
#include "Error.h"
#include "Init.h"
#include "ReadingArguments.h"

void ReadArgs(int argc, char *argv[], init *I)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-d")) {
			int l;
			sscanf(argv[++i], "%d", &l);
			I->SetDebugLevel(l);
		} else if (!strcmp(argv[i], "-m")) {
			I->SetInputFileName(argv[++i]);
		} else
			cout << INPUT_DATA_ERR;//DebugP(0, INPUT_DATA_ERR << argv[i]);
	}
}
