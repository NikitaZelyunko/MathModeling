#ifndef ABSTRACT_SOLVER_H
#define ABSTRACT_SOLVER_H

#include <iostream>

/* R - result type, S - scalar type (epsilon type) */
template<class R, class S>
class AbstractSolver {
private:
    R result;
protected:

    bool isInit = false;
    S epsilon;

    virtual void solverNotInited() const {
        std::cout<<std::endl<<"Solver not inited"<<std::endl;
    }

    virtual void solverNotSolved() const {
        std::cout<<std::endl<<"Solver not solved. Run solve() method before ask result."<<std::endl;
    }

    virtual void resultIsNotInited() const {
        std::cout<<std::endl<<"result is not inited"<<std::endl;
    }

    /* 
        Define this in your class (remember about function signature)
        const R computeResult() const;
    */
    virtual const R computeResult() const {return this->getResult();}

    const R setResult(const R& newValue) {
        isInit = true;
        result = newValue;
        return result;
    }

    AbstractSolver() {
        this->epsilon = 0;
    }

    AbstractSolver(S epsilon) {
        this->epsilon = epsilon;
    }

public:
    const R solve() {
        this->setResult(computeResult());
        return this->getResult();
    };
    
    ~AbstractSolver(){};

    bool isSolved() const {
        return this->isInit;
    }

    const R getResult() const {
        if(!isInit) {
           resultIsNotInited();
        }
        return this->result;
    }

    const S getEpsilon() const {
        return this->epsilon;
    }

};

#endif