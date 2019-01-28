#pragma once

template <class T>
class Singleton
{
	static T* _self;
	static int _refcount;
protected:
	Singleton(){}
	virtual ~Singleton(){_self = NULL;}
public:
	static T* Instance();
	void FreeInst();
};

template <class T>
T*  Singleton<T>::_self = NULL;

template <class T>
int  Singleton<T>::_refcount = 0;
