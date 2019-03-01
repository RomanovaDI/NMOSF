#pragma once

#define FUNC_START "Function started"
#define FUNC_FINISH "Function finished"
#define MEM_ERR "Memory error:"
#define INIT_DATA_ERR "Initial data error:"
#define FILE_OPEN_ERR "No such file:"
#define FILE_ERR "Error file:"
#define FILE_DATA_ERR "Error data in file:"
#define INPUT_DATA_ERR "Error input data:"
#define MAP_FILE_ERR "Error file of map:"
#define VECT_COMP_ERR "Vector component error: No component"
#define INVALID_CHAR_ERR "Invalid simbol"

#define DEBUGLVL 10
#define EXITLVL 0

#define DEBUGP(level, args...) DebugP(level, __FILE__, ":", __LINE__, ":",  __func__, ":", args)

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
void DebugP(int level, Args ... args)
{
	if (level <= DEBUGLVL) {
		do_it(args ...);
		cout << endl;
		if (level <= EXITLVL)
			exit(0);
	}
}
