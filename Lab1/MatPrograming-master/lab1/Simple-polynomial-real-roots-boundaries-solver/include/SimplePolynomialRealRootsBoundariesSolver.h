#ifndef SIMPLE_POLYNOMIAL_REAL_ROOTS_BOUNDARIES_SOLVER_H
#define SIMPLE_POLYNOMIAL_REAL_ROOTS_BOUNDARIES_SOLVER_H

#include "AbstractSolver.h"
#include "Point.h"
#include "SimplePolynomial.h"

const short int SIMPLE_METHOD = 0;
const short int NEWTON_METHOD = 1;

template <class S>
class SimplePolynomialRealRootsBoundariesSolver : public AbstractSolver<Point<Point<S>>, S> {
    typedef Point<Point<S>> BoundariesResult;
protected:
    Point<S> polinomialKoeffs;
    short int methodType;

    const BoundariesResult simpleMethod() const {
        S* polynomialKoeffs = Point<S>::toArray(this->polinomialKoeffs);
        S positiveCheck = findMin(polynomialKoeffs, this->polinomialKoeffs.length());

        this->polinomialKoeffs.forEach([&](int i, const Point<S>& point) -> void {
            polynomialKoeffs[i] = abs(polynomialKoeffs[i]);
        });

        S A = findMax(polynomialKoeffs, this->polinomialKoeffs.length() - 1);
        S range = 1 + A / abs(this->polinomialKoeffs[this->polinomialKoeffs.length()-1]);

        Point<S> positiveBoundaries;
        if(positiveCheck < 0) {
            positiveBoundaries = Point<S>(2);
            positiveBoundaries[0] = 0;
            positiveBoundaries[1] = range;
        }
        Point<S> negativeBoundaries(2);
        negativeBoundaries[0] = -range;
        negativeBoundaries[1] = 0;

        Point<Point<S>> result(2);
        result[0] = negativeBoundaries;
        result[1] = positiveBoundaries;
        return result;
    }

    /* TODO move part of this to single solver **/
    const S secantMethod(const SimplePolynomial<S>& polynomial, const S start, const S step) const {
        S x1, x2, res, y, y1, y2, h = step;
        x1 = start;
        x2 = start + h;
        y1 = polynomial(x1);
        y2 = polynomial(x2);
        S forecast = x1;
        S forecastValue = y1;

        for(; forecastValue < -this->getEpsilon();) {
            if(abs(y1 - y2) > this->getEpsilon()) {
                forecast = (x2 * y1 - x1*y2)/(y1 - y2);
            }
            forecastValue = polynomial(forecast);
            if(forecastValue < -this->getEpsilon()) {
                x1 = forecast;
                x2 = x1 + step;
                y1 = polynomial(x1);
                y2 = polynomial(x2);
            }
        }
        x2 = forecast;
        /* TODO this **/
        if(forecastValue > this->getEpsilon()) {
            if(abs(x2-x1)> this->getEpsilon()){
                do
                {
                    y = res;
                    res = x2 - ((x2 - x1) / (polynomial(x2) - polynomial(x1))) * polynomial(x2);
                    x1 = x2;
                    x2 = res;
                }
                while (abs(y - res) >= this->getEpsilon());
            } else {
                res = x2;
            }
        } else {
            res = x2;
        }
        /* TODO this **/
        return res;
    }

    const S newtonMethodStep(const SimplePolynomial<S>& polynomial) const {
        SimplePolynomial<S> bufPolynomial(polynomial);
        Point<SimplePolynomial<S>> derivatives(bufPolynomial.getDegree() + 1);
        derivatives.forEach([&](int i, const Point<SimplePolynomial<S>>& point) -> void {
            if(i == 0) {
                point[0] = bufPolynomial;
            } else {
                point[i] = bufPolynomial.derivative();
            }
        });
        return derivatives.template reduceReverse<S>(0.0, [&](const S acc, int i, const Point<SimplePolynomial<S>> point) -> const S {
           S step = 1;
           SimplePolynomial<S> *f = &point[i];
           S res = secantMethod(*f, acc, step);
           return res;
        });
    }

    const BoundariesResult newtonMethod() const {
        S first = newtonMethodStep(SimplePolynomial<S>(polinomialKoeffs));
        std::cout<<"FIRST: "<<first<<std::endl;
        return Point<Point<S>>(2, Point<S>(2, 3));
    }

    const BoundariesResult computeResult() const {
        switch (methodType) {
            case SIMPLE_METHOD: {
                return simpleMethod();
            }
            case NEWTON_METHOD: {
                return newtonMethod();
            }
        }
        return Point<Point<S>>(2, Point<S>(2, 3));
    }

public:
    using AbstractSolver<BoundariesResult, S>::AbstractSolver;

    SimplePolynomialRealRootsBoundariesSolver(int n, S* koeffs, const short int methodType, const S epsilon):AbstractSolver<BoundariesResult, S>(epsilon) {
        this->methodType = methodType;
        polinomialKoeffs = Point<S>(n, koeffs);
    }

    SimplePolynomialRealRootsBoundariesSolver(const SimplePolynomial<S>& polynomial, const short int methodType, const S epsilon):AbstractSolver<BoundariesResult, S>(epsilon) {
        this->methodType = methodType;
        polinomialKoeffs = polynomial.getKoeffs();
    }
};

void TestSimplePolynomialRealRootsBoundariesSolver();

#endif //SIMPLE_POLYNOMIAL_REAL_ROOTS_BOUNDARIES_SOLVER_H
