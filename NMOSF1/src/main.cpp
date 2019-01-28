#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
//#include "MathModel.h"
//nclude "Error.h"
//#include "Init.h"
//#include "Singleton.h"
//#include "Mesh.h"

template <class T>
class Singleton
{
	static T* _self;
	static int _refcount;
protected:
	Singleton(){}
	virtual ~Singleton(){_self = NULL;}
public:
	static T* Instance();
	void FreeInst();
};

template <class T>
T*  Singleton<T>::_self = NULL;

template <class T>
int  Singleton<T>::_refcount = 0;

template <class T>
T*  Singleton<T>::Instance()
{
	if(!_self)
		_self=new T;
	_refcount++;
	return _self;
}

template <class T>
void  Singleton<T>::FreeInst()
{
	if(--_refcount==0)
		delete this;
}

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

init::init() {}

void init::ReadInputFile(char name[100])
{
	strcpy(InputFileName, name);
	cout << InputFileName << endl;
	cout << "blaaaaaaaaaaaaaa" << endl;
	FILE *InputFile = fopen(InputFileName, "r");
	if (InputFile == NULL)
		cout << "Error input file" << endl;//DebugCout(0, FILE_OPEN_ERR + ": " + InputFileName);
	cout << "File was opened" << endl;
	char str[300];
	char key[100];
	while (fgets(str, 300, InputFile)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "DEBUG")) {
				int debugLevel;
				if (! sscanf(str, "%s %d", key, &debugLevel))
					cout << "error debug level" << endl;//DebugCout(0, FILE_DATA_ERR + ": error debug level in file " + InputFileName);
				valDebugLevel = debugLevel;
			} else if (! strcmp(key, "MapName")) {
				char mapName[100];
				if (! sscanf(str, "%s %s", key, mapName))
					cout << "error mapName" << endl;//DebugCout(0, FILE_DATA_ERR + ": error MapName in file " + InputFileName);
				strcpy(MapName, mapName);
			} else if (! strcmp(key, "RegionFileName")) {
				char regionFileName[100];
				if (! sscanf(str, "%s %s", key, regionFileName))
					cout << "error regionFileName" << endl;//DebugCout(0, FILE_DATA_ERR + ": error RegionFileName in file " + InputFileName);
				strcpy(RegionFileName, regionFileName);
			} else
				cout << "unknown tag" << endl;//DebugCout(0, FILE_DATA_ERR + ": unknown tag " + key + " in file " + InputFileName);	
		} else
			cout << "error key" << endl;
	}
	fclose(InputFile);
}

void init::PrintInfo()
{
	cout << InputFileName << endl;
	cout << MapName << endl;
	cout << RegionFileName << endl;
	cout << valDebugLevel << endl;
}

int init::DebugLevel()
{
	return valDebugLevel;
}

int main(int argc, char* argv[])
{
	cout << "starting" << endl;
	if (argc != 2)
		cout << "not enouth arguments" << endl;//DebugCout(0, INIT_DATA_ERR);
	cout << "Correct numbeg of arguments" << endl;
	init *Init = init::Instance();
	cout << "initializing init" << endl;
	Init->ReadInputFile(argv[1]);
	cout << "input file was read" << endl;
	Init->PrintInfo();
	//mesh myMesh(initData);
    return 0;
}
