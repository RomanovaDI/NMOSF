#pragma once
#include <stdio.h>

#define MEM_ERR "Memory error"
#define INIT_DATA_ERR "Initial data error"
#define FILE_ERR "Error file"
#define FILE_DATA_ERR "Error data in file"

#define DEBUG 10
#define DebugCout(level, args)\
	{if (level <= DEBUG)\
		cout << __FILE__ << ":" << __func__ << ":" << __LINE__ << ": " << args << endl}
