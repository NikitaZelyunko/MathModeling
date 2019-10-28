#ifndef FUNCTION_H
#define FUNCTION_H

#include "Point.h"

using namespace std;

template<class T>
class Function
{
public:
    virtual const T operator() (const Point<T>& p) { return p[0]; };
    Function() = default;
    virtual  ~Function() = default;
};

#endif