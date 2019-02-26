#pragma once

class variable
{
private:
	char varName[50];
	double *varVal;
	int varVect;
public:
	variable() {}
	~variable() {}
	void setName(char[50]);
	void getName(char *);
	void setVect(int);
	int getVect();
};

class mathModel
{
private:
	variable *varList;
	int varListLen;
	void readVarList(char[100]);
	void readVarEqn(char[100]);
public:
	mathModel() {}
	void readMathModel(char[100]);
	~mathModel();
};
