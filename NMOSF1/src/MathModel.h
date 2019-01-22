#pragma once

class variable
{
private:
	char varName[50];
	double *varVal;
public:
	variable(int);
	void setName(char[50]);
	char* Name();
	double Value(int, int, int);
	double* Address(int, int, int);
}
