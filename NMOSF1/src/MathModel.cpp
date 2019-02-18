#include <iostream>
#include <string.h>
using namespace std;
#include "Error.h"
#include "MathModel.h"
 
void variable::setName(char name[50])
{
	strcpy(varName, name);
}

void variable::getName(char *name)
{
	strncpy(name, varName, 50);
}

void mathModel::readMathModel(char fileName[100])
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[500];
	char varList_tmp[100][50];
	varListLen = 0;
	int flag = 0;
	char key[50];
	while (fgets(str, 500, f)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "StartVarList"))
				flag++;
			else if (! strcmp(key, "FinishVarList"))
				break;
			else if (flag == 1) {
				strcpy(varList_tmp[varListLen], key);
				varListLen++;
				if (varListLen == 100)
					DEBUGP(0, FILE_DATA_ERR, "Too many variables. Check \"FinishVarList\" tag.");
			} else
				DEBUGP(0, FILE_DATA_ERR, "no tag \"StartVarList\" in file", fileName);
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
	}
	varList = new variable[varListLen];
	for (int i = 0; i < varListLen; i++)
		varList[i].setName(varList_tmp[i]);
	DEBUGP(1, FUNC_FINISH);
}

mathModel::~mathModel()
{
	delete[] varList;
}
