#pragma once

#define MEM_ERR "Memory error :"
#define INIT_DATA_ERR "Initial data error :"
#define FILE_OPEN_ERR "No such file :"
#define FILE_ERR "Error file :"
#define FILE_DATA_ERR "Error data in file :"
#define INPUT_DATA_ERR "Error input data :"

class error
{
protected:
	int valDebugLevel;
	template<typename Arg>
	void do_print(Arg a);
	template<typename Arg>
	void do_print<Arg>(Arg a);
	template<typename Arg>
	void do_it(Arg a);
	template<typename Arg, typename... Args>
	void do_it(Arg arg, Args... args);
public:
	error();
	void SetDebugLevel(int);
	int GetDebugLevel();
	template<typename... Args>
	void DebugP(int, const string, Args... args);
};
