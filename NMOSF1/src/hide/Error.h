#pragma once

#define MEM_ERR "Memory error"
#define INIT_DATA_ERR "Initial data error"
#define FILE_OPEN_ERR "No such file"
#define FILE_ERR "Error file"
#define FILE_DATA_ERR "Error data in file"

//#define DebugCout(level, args)\
//	if (level <= Init.DebugLevel()) {\
//		cout << __FILE__ << ":" << __func__ << ":" << __LINE__ << ": " << args << endl;\
//		exit(0);}

class error
{
private:
	int valDebugLevel;
public:
	error(init);
	DebugP(int, string);
};
