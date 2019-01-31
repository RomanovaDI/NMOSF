#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
#include "Error.h"

std::stringstream MSG;

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

template<typename Arg>
void error::do_print(Arg a)
{
    std::cout << " " << a << " ";
}

template<typename Arg>
void error::do_print<Arg>(Arg a)
{
    std::cout << " " << a << " ";
}

template<typename Arg>
void error::do_it(Arg a)
{
    do_print(a);
}

template<typename Arg, typename... Args>
void error::do_it(Arg arg, Args... args)
{
    do_print(arg);
	do_it(args...);
}

template<typename... Args>
void error::DebugP(int level, const string err, Args... args)
{
	if (level <= valDebugLevel) {
		cout << err;
		do_it(args...);
		exit(0);
	}
}
