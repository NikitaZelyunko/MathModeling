#ifndef ONE_DIMENSIONAL_FUNCTION_H
#define ONE_DIMENSIONAL_FUNCTION_H
#include "Function.h"

template <class T>
class OneDimensionalFunction : public Function<T>
{
public:
    using Function<T>::Function;
    virtual const T operator() (const T& point) { return point; };
    virtual ~OneDimensionalFunction();
};

template<class T>
OneDimensionalFunction<T>::~OneDimensionalFunction() {

}

#endif