#include <iostream>
using namespace std;
#include "MathModel.h"
#include "Error.h"
 
variable::variable(int MeshSize)
{
	if ((varVal = (double *) malloc(sizeof(double) * MeshSize)) == NULL)
		DebugCout(0, MEM_ERR);
	return;
}

void variable::setName(char name[50])
{
	varName = name;
	return;
}

char* variable::Name()
{
	return varName;
}

double variable::Value(int i)
{
	return varVal[i];
}
 
void CppStudio::message() // функция (метод класса) выводящая сообщение на экран
{
 cout << "nwebsite: cppstudio.comntheme: Classes and Objects in C + +n";
}
 
void CppStudio::setDate(int date_day, int date_month, int date_year) // установка даты в формате дд.мм.гг
{
 day   = date_day; // инициализация день
 month = date_month; // инициализация месяц
 year  = date_year; // инициализация год
}
 
void CppStudio::getDate() // отобразить текущую дату
{
 cout << "date: " << day << "." << month << "." << year << endl;
}
