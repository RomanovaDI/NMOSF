#pragma once

#define MEM_ERR "Memory error"
#define INIT_DATA_ERR "Initial data error"
#define FILE_OPEN_ERR "No such file"
#define FILE_ERR "Error file"
#define FILE_DATA_ERR "Error data in file"
#define INPUT_DATA_ERR "Error input data"

class error
{
protected:
	int valDebugLevel;
public:
	error();
	void SetDebugLevel(int);
	int GetDebugLevel();
	void DebugP(int, string);
};
