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

void variable::setVect(int vect)
{
	varVect = vect;
}

int variable::getVect()
{
	return varVect;
}

void equation::setTag(char tag[50])
{
	strcpy(eqnTag, tag);
}

void equation::getTag(char *tag)
{
	strcpy(tag, eqnTag);
}

void equation::setRecord(char record[500])
{
	strcpy(eqnRecord, record);
}

void equation::getRecord(char *record)
{
	strcpy(record, eqnRecord);
}

void mathModel::readVarList(char fileName[100])
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[500];
	char varList_tmp[100][50];
	int vect[100];
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
				if ((strlen(key) > 2) && (key[strlen(key) - 2] == '[') && (key[strlen(key) - 1] == ']')) {
					strncpy(varList_tmp[varListLen], key, strlen(key) - 2);
					vect[varListLen++] = 1;
				} else {
					strcpy(varList_tmp[varListLen], key);
					vect[varListLen++] = 0;
				}
				if (varListLen == 100)
					DEBUGP(0, FILE_DATA_ERR, "Too many variables. Check \"FinishVarList\" tag.");
			} else
				continue;
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
	}
	varList = new variable[varListLen];
	for (int i = 0; i < varListLen; i++) {
		varList[i].setName(varList_tmp[i]);
		varList[i].setVect(vect[i]);
	}
	calcVarListLenVect();
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::calcVarListLenVect()
{
	DEBUGP(1, FUNC_START);
	int tmp;
	//for (int i = 0, varListLenVect = 0; i < varListLen; i++)
	for (int i = 0; i < varListLen; i++)
		if (tmp = varList[i].getVect())
			eqnListLen += dimension;
		else
			eqnListLen++;
	printf("eqnListLen = %d\n", eqnListLen);
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::readVarEqn(char fileName[100])
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[500];
	char eqnListTag_tmp[100][50];
	char eqnList_tmp[100][500];
	int listLen_tmp = 0;
	int flag = 0;
	char key[50];
	while (fgets(str, 500, f)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			printf("%s\n", key);
			if (! strcmp(key, "StartVarEqn"))
				flag++;
			else if (! strcmp(key, "FinishVarEqn"))
				break;
			else if (flag == 1) {
				char *colon_pointer = strchr(str, ':');
				if (colon_pointer == NULL)
					DEBUGP(0, FILE_DATA_ERR, "No expression \"VarName:\".");
				int colon = colon_pointer - str;
				char *square_br_pointer = strchr(str, '[');
				int square_br = square_br_pointer - str;
				int var_len;
				if (square_br < colon)
					var_len = square_br;
				else
					var_len = colon;
				strncpy(eqnListTag_tmp[listLen_tmp], key, var_len);
				strcpy(eqnList_tmp[listLen_tmp], str);
				listLen_tmp++;
				if (listLen_tmp > eqnListLen)
					DEBUGP(0, FILE_DATA_ERR, "Too many equations. Check \"FinishVarList\" tag.");
			} else
				continue;
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
	}
	eqnList = new equation[eqnListLen];
	for (int i = 0; i < eqnListLen; i++) {
		eqnList[i].setTag(eqnListTag_tmp[i]);
		eqnList[i].setRecord(eqnList_tmp[i]);
	}
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::readMathModel(char fileName[100], int dim)
{
	dimension = dim;
	DEBUGP(1, FUNC_START);
	readVarList(fileName);
	readVarEqn(fileName);
	for (int i = 0; i < varListLen; i++) {
		char var_tmp[50];
		varList[i].getName(var_tmp);
		printf("%s[%d]\n", var_tmp, varList[i].getVect() * dimension);
	}
	for (int i = 0; i < eqnListLen; i++) {
		char tag_tmp[50];
		char rec_tmp[500];
		eqnList[i].getTag(tag_tmp);
		eqnList[i].getRecord(rec_tmp);
		printf("%s|||%s\n", tag_tmp, rec_tmp);
	}
	DEBUGP(1, FUNC_FINISH);
}

mathModel::~mathModel()
{
	delete[] varList;
	delete[] eqnList;
}
