#pragma once
#include <stdio.h>

#define MEM_ERR "Memory error"

#define DEBUG 10
#define DebugCout(level, args)\
	if (level <= DEBUG)\
		cout << __FILE__ << ":" << __func__ << ":" << __LINE__ << ": " << args << endl 
