#pragma once

#define TAGLEN 50
#define LINELEN 500
#define MAXLIST 100

int deleteSpLF(char* str, int len);
int whichSimbol(char ch);

class variable
{
private:
	char varName[TAGLEN];
	double *varVal;
	int varVect;
public:
	variable();
	~variable() {}
	void setName(char[TAGLEN]);
	void getName(char*);
	void setVect(int);
	int getVect();
};

class equation
{
private:
	char eqnTag[TAGLEN];
	int eqnTagComponent;
	char eqnRecord[LINELEN];
public:
	equation() {}
	~equation() {}
	void setTag(char[TAGLEN]);
	void getTag(char*);
	void setRecord(char[LINELEN]);
	void getRecord(char*);
	void setTagComponent(int);
	int getTagComponent();
};

class constant
{
private:
	char constName[TAGLEN];
	int constDimension = 0;
	double constValue[3] = {0, 0, 0};
public:
	constant() {}
	~constant() {}
	void setName(char[TAGLEN]);
	void getName(char*);
	void setValue(double);
	void setValue(int, double);
	double getValue();
	double getValue(int);
};

class mathModel
{
private:
	int dimension;
	variable *varList;
	equation *eqnList;
	constant *constList;
	int varListLen;
	int eqnListLen;
	int constListLen;
	void calcEqnListLen();
	void readVarList(char*);
	void readVarEqnList(char*);
	void readConstList(char*);
public:
	mathModel() {}
	void readMathModel(char*, int);
	~mathModel();
};
