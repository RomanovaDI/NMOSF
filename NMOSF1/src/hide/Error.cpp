#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
#include "Error.h"

error::error()
{
	valDebugLevel = 10;
}

void error::SetDebugLevel(int l)
{
	valDebugLevel = l;
}

int error::GetDebugLevel()
{
	return valDebugLevel;
}
