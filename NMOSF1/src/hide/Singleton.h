#pragma once

template <class T>
class Singleton
{
	static T* _self;
protected:
	Singleton(){}
	virtual ~Singleton(){}
public:
	static T* Instance();
	void FreeInst();
};
