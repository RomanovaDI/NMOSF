#pragma once

#define MEM_ERR "Memory error :"
#define INIT_DATA_ERR "Initial data error :"
#define FILE_OPEN_ERR "No such file :"
#define FILE_ERR "Error file :"
#define FILE_DATA_ERR "Error data in file :"
#define INPUT_DATA_ERR "Error input data :"
#define MAP_FILE_ERR "Error file of map :"

#define DEBUGLVL 10

template<typename Arg>
void do_print(Arg a)
{
 std::cout << a << " ";
}

template<typename Arg>
void do_it(Arg a)
{
	do_print(a);
}

template<typename Arg, typename ... Args>
void do_it(Arg arg, Args ... args)
{
	do_print(arg);
	do_it(args ...);
}

template<typename ... Args>
void DEBUGP(int level, Args ... args)
{
	if (level <= DEBUGLVL) {
		do_it(args ...);
		cout << endl;
		exit(0);
	}
}
