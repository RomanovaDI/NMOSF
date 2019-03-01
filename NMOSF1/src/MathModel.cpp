#include <iostream>
#include <string.h>
using namespace std;
#include "Error.h"
#include "MathModel.h"

variable::variable()
{
	memset(varName, 0, TAGLEN * sizeof(char));
}

void variable::setName(char name[TAGLEN])
{
	strcpy(varName, name);
}

void variable::getName(char *name)
{
	strncpy(name, varName, TAGLEN);
}

void variable::setVect(int vect)
{
	varVect = vect;
}

int variable::getVect()
{
	return varVect;
}

void equation::setTag(char tag[TAGLEN])
{
	strcpy(eqnTag, tag);
}

void equation::getTag(char *tag)
{
	strcpy(tag, eqnTag);
}

void equation::setRecord(char record[LINELEN])
{
	strcpy(eqnRecord, record);
}

void equation::getRecord(char *record)
{
	strcpy(record, eqnRecord);
}

void equation::setTagComponent(int comp)
{
	eqnTagComponent = comp;
}

int equation::getTagComponent()
{
	return eqnTagComponent;
}

int deleteSpLF(char* str, int len)
{
	for (int i = len - 1; i >= 0; i--)
		if (str[i] == '\n') {
			str[i] = '\0';
			len--;
		}
	for (int i = 0; i < len; i++)
		if (str[i] == ' ') {
			for (int j = i; j < len - 1; j++)
				str[j] = str[j + 1];
			str[len - 1] = '\0';
			len--;
		}
}

int whichSimbol(char ch)
{
	int tmp = ch;
	if ((tmp >= '0') && (tmp <= '9'))
		return 1;
	if (((tmp >= 'A') && (tmp <= 'Z')) || ((tmp >= 'a') && (tmp <= 'z')))
		return 2;
	if (tmp == '_')
		return 3;
	if (tmp == '=')
		return 4;
	if ((tmp == '*') || (tmp == '/') || (tmp == '+') || (tmp == '-') || (tmp == '^') || (tmp == '?') || tmp == ':')
		return 5;
	if (tmp == '.')
		return 6;
	if (tmp == '(')
		return 7;
	if (tmp == ')')
		return 8;
	if (tmp == '[')
		return 9;
	if (tmp == ']')
		return 10;
	return 0;
}

void constant::setName(char name[TAGLEN])
{
	strcpy(constName, name);
}

void constant::getName(char *name)
{
	strcpy(name, constName);
}

void constant::setValue(double value)
{
	constDimension = 1;
	constValue[0] = value;
}

void constant::setValue(int comp, double value)
{
	if (comp >= constDimension)
		constDimension = comp + 1;
	constValue[comp] = value;
}

double constant::getValue()
{
	return constValue[0];
}

double constant::getValue(int comp_tmp)
{
	if (comp_tmp >= 3)
		DEBUGP(0, VECT_COMP_ERR, comp_tmp);
	return constValue[comp_tmp];
}

void mathModel::readVarList(char *fileName)
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[LINELEN];
	memset(str, 0, LINELEN * sizeof(char));
	char varList_tmp[MAXLIST][TAGLEN];
	int vect[MAXLIST];
	varListLen = 0;
	int flag = 0;
	char key[TAGLEN];
	while (fgets(str, LINELEN, f)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "StartVarList"))
				flag++;
			else if (! strcmp(key, "FinishVarList"))
				break;
			else if (flag == 1) {
				deleteSpLF(str, strlen(str));
				if ((strlen(str) > 2) && (str[strlen(str) - 2] == '[') && (str[strlen(str) - 1] == ']')) {
					strncpy(varList_tmp[varListLen], str, strlen(str) - 2);
					vect[varListLen++] = 1;
				} else {
					strcpy(varList_tmp[varListLen], str);
					vect[varListLen++] = 0;
				}
				if (varListLen == MAXLIST)
					DEBUGP(0, FILE_DATA_ERR, "Too many variables. Check \"FinishVarList\" tag.");
			} else
				continue;
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
		memset(str, 0, LINELEN * sizeof(char));
	}
	varList = new variable[varListLen];
	for (int i = 0; i < varListLen; i++) {
		for (int j = 0; j < strlen(varList_tmp[i]); j++)
			if (!whichSimbol(varList_tmp[i][j]))
				DEBUGP(0, INVALID_CHAR_ERR, varList_tmp[i][j]);
		varList[i].setName(varList_tmp[i]);
		varList[i].setVect(vect[i]);
	}
	calcEqnListLen();
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::calcEqnListLen()
{
	DEBUGP(1, FUNC_START);
	int tmp;
	for (int i = 0; i < varListLen; i++)
		if (tmp = varList[i].getVect())
			eqnListLen += dimension;
		else
			eqnListLen++;
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::readVarEqnList(char *fileName)
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[LINELEN];
	char eqnListTag_tmp[MAXLIST][TAGLEN];
	char eqnList_tmp[MAXLIST][LINELEN];
	int eqnListTagComp_tmp[MAXLIST];
	eqnListLen = 0;
	int flag = 0;
	char key[TAGLEN];
	while (fgets(str, LINELEN, f)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "StartVarEqnList"))
				flag++;
			else if (! strcmp(key, "FinishVarEqnList"))
				break;
			else if (flag == 1) {
				deleteSpLF(str, strlen(str));
				char *colon_pointer = strchr(str, ':');
				if (colon_pointer == NULL)
					DEBUGP(0, FILE_DATA_ERR, "No expression \"VarName:\".");
				int colon = colon_pointer - str;
				char *square_br_pointer = strchr(str, '[');
				int square_br = square_br_pointer - str;
				int var_len;
				if (square_br < colon) {
					var_len = square_br;
					if (square_br != colon - 3)
						DEBUGP(0, FILE_DATA_ERR, "Error tag of eqation");
					if (str[colon - 1] != ']')
						DEBUGP(0, FILE_DATA_ERR, "Error tag of eqation");
					int comp_tmp = str[square_br + 1] - '0';
					if ((comp_tmp < 0) || (comp_tmp >= dimension))
						DEBUGP(0, FILE_DATA_ERR, "Error component", comp_tmp);
					eqnListTagComp_tmp[eqnListLen] = comp_tmp;
				} else
					var_len = colon;
				strncpy(eqnListTag_tmp[eqnListLen], key, var_len);
				strcpy(eqnList_tmp[eqnListLen++], &str[colon + 1]);
				if (eqnListLen > eqnListLen)
					DEBUGP(0, FILE_DATA_ERR, "Too many equations. Check \"FinishVarList\" tag.");
			} else
				continue;
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
	}
	printf("1\n");
	printf("%d\n", eqnListLen);
	eqnList = new equation[eqnListLen];
	for (int i = 0; i < eqnListLen; i++) {
		eqnList[i].setTag(eqnListTag_tmp[i]);
		eqnList[i].setRecord(eqnList_tmp[i]);
		eqnList[i].setTagComponent(eqnListTagComp_tmp[i]);
	}
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::readConstList(char *fileName)
{
	DEBUGP(1, FUNC_START);
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
		DEBUGP(0, FILE_OPEN_ERR, fileName);
	char str[LINELEN];
	char constList_tmp[MAXLIST][TAGLEN];
	double constListValue_tmp[MAXLIST];
	memset(constListValue_tmp, 0, sizeof(double) * MAXLIST);
	int constListComp_tmp[MAXLIST];
	memset(constListComp_tmp, 0, sizeof(int) * MAXLIST);
	int constListLen_tmp = 0;
	int flag = 0;
	char key[TAGLEN];
	while (fgets(str, LINELEN, f)) {
		if (str[0] == '#')
			continue;
		else if (sscanf(str, "%s", key)) {
			if (! strcmp(key, "StartConstList"))
				flag++;
			else if (! strcmp(key, "FinishConstList"))
				break;
			else if (flag == 1) {
				deleteSpLF(str, strlen(str));
				char *equality_pointer = strchr(str, '=');
				if (equality_pointer == NULL)
					DEBUGP(0, FILE_DATA_ERR, "No expression \"const_name = value\".");
				int equality = equality_pointer - str;
				int key_end = equality;
				if ((str[equality - 1] == ']') && (str[equality - 3] == '[')) {
					if (whichSimbol(str[equality - 2]) != 1)
						DEBUGP(0, INVALID_CHAR_ERR, str[equality - 2]);
					else
						constListComp_tmp[constListLen_tmp] = str[equality - 2] - '0';
					key_end = equality - 3;
				}
				strncpy(constList_tmp[constListLen_tmp], str, key_end);
				constListValue_tmp[constListLen_tmp++] = strtod(&str[equality + 1], NULL);
				if (constListLen_tmp == MAXLIST)
					DEBUGP(0, FILE_DATA_ERR, "Too many constants. Check \"FinishConstList\" tag.");
			} else
				continue;
		} else
			DEBUGP(0, FILE_DATA_ERR, "in file", fileName, "in line:", str);
	}
	constListLen = 0;
	for (int i = 0; i < constListLen_tmp; i++) {
		int f = 0;
		for (int j = i + 1; j < constListLen_tmp; j++)
			if (strcmp(constList_tmp[i], constList_tmp[j]) == 0) {
				f++;
				break;
			}
		if (!f)
			constListLen++;
	}
	constList = new constant[constListLen];
	for (int i = 0; i < constListLen_tmp; i++) {
		int j;
		for (j = 0; j < i; j++)
			if (strcmp(constList_tmp[i], constList_tmp[j]) == 0)
				break;
		constList[j].setName(constList_tmp[i]);
		constList[j].setValue(constListComp_tmp[i], constListValue_tmp[i]);
	}
	DEBUGP(1, FUNC_FINISH);
}

void mathModel::readMathModel(char *fileName, int dim)
{
	dimension = dim;
	DEBUGP(1, FUNC_START);
	readVarList(fileName);
	readVarEqnList(fileName);
	readConstList(fileName);
	for (int i = 0; i < varListLen; i++) {
		char var_tmp[TAGLEN];
		varList[i].getName(var_tmp);
		printf("%s[%d]\n", var_tmp, varList[i].getVect() * dimension);
	}
	for (int i = 0; i < eqnListLen; i++) {
		char tag_tmp[TAGLEN];
		char rec_tmp[LINELEN];
		eqnList[i].getTag(tag_tmp);
		eqnList[i].getRecord(rec_tmp);
		int comp_tmp = eqnList[i].getTagComponent();
		printf("%s[%d]: %s\n", tag_tmp, comp_tmp, rec_tmp);
	}
	for (int i = 0; i < constListLen; i++) {
		char name_tmp[TAGLEN];
		constList[i].getName(name_tmp);
		double value_tmp[3];
		value_tmp[0] = constList[i].getValue(0);
		value_tmp[1] = constList[i].getValue(1);
		value_tmp[2] = constList[i].getValue(2);
		printf("%s = {%lf, %lf, %lf}\n", name_tmp, value_tmp[0], value_tmp[1], value_tmp[2]);
	}
	DEBUGP(1, FUNC_FINISH);
}

mathModel::~mathModel()
{
	delete[] varList;
	delete[] eqnList;
}
