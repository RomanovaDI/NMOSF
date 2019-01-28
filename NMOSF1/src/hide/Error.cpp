#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "Error.h"
#include "Init.h"

error::error(init Init)
{
	valDebugLevel = Init.valDebugLevel;
}

error::DebugP(int level, string message)
{
	if (level <= valDebugLevel) {
		cout << message;
		exit(0);
	}
}
