#include <iostream>
#include <fstream>
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

void error::DebugP(int level, string message)
{
	if (level <= valDebugLevel) {
		cout << message;
		exit(0);
	}
}
