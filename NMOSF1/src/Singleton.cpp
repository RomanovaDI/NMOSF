#include <iostream>
#include <stdio.h>
using namespace std;
#include "Singleton.h"

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
