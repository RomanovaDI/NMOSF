#include <iostream>
#include <string>
#include <string.h>
using namespace std;
#include "ReadingArguments.h"
#include "Init.h"
#include "Error.h"

int ReadArgs(int argc, char *argv[], init *I)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-dl")) {
			int l;
			sscanf(argv[++i], "%d", l);
			I->SetDebugLevel(l);
		} else if (!strcmp(argv[i], "-m")) {
			I->SetMapName(argv[++i]);
		} else
			DebugP(0, INPUT_DATA_ERR + argv[i]);
	}
}
