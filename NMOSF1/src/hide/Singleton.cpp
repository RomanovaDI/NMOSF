#include <iostream>
using namespace std;
#include "Singleton.h"

template <class T>
T*  Singleton<T>::_self = NULL;

template <class T>
T*  Singleton<T>::Instance()
{
	if(!_self)
		_self = new T;
	return _self;
}

template <class T>
void  Singleton<T>::FreeInst()
{
	delete this;
}
