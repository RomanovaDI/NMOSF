#pragma once

#include "Singleton.h"

//class error;

class init : public Singleton<init>
{
private:
	char InputFileName[100];
	char MapName[100];
	char RegionFileName[100];
	int valDebugLevel;
public:
	void ReadInputFile(char [100]);
	void PrintInfo();
	int DebugLevel();
//	friend error;
protected:
	init();
	friend class Singleton<init>;
};
