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
	void getName(char*);
	void setVect(int);
	int getVect();
};

class equation
{
private:
	char eqnTag[50];
	char eqnRecord[500];
public:
	equation() {}
	~equation() {}
	void setTag(char[50]);
	void getTag(char*);
	void setRecord(char[500]);
	void getRecord(char*);
};

class mathModel
{
private:
	int dimension;
	variable *varList;
	equation *eqnList;
	int varListLen;
	int eqnListLen = 0;
	void calcVarListLenVect();
	void readVarList(char[100]);
	void readVarEqn(char[100]);
public:
	mathModel() {}
	void readMathModel(char[100], int);
	~mathModel();
};
