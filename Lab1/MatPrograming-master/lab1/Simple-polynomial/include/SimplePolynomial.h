#ifndef SIMPLE_POLYNOMIAL_H
#define SIMPLE_POLYNOMIAL_H

#include "Point.h"
#include "OneDimensionalFunction.h"
#include "Utils.h"

template <class T>
class SimplePolynomial: public OneDimensionalFunction<T> {
private:
    Point<T> koeffs;

public:
    using OneDimensionalFunction<T>::OneDimensionalFunction;

    SimplePolynomial(int n):SimplePolynomial<T>() {
        this->koeffs = Point<T>(n);
    }

    SimplePolynomial(const Point<T>& koeffs): SimplePolynomial<T>() {
        this->koeffs = Point<T>(koeffs);
    }
    SimplePolynomial(const SimplePolynomial<T>& polinom): SimplePolynomial<T>(polinom.koeffs) {}

    SimplePolynomial(int n, T* koeffs) {
        this->koeffs = Point<T>(n, koeffs);
    }
    SimplePolynomial(int n, const T& filler) {
        this->koeffs = Point<T>(n, filler);
    }

    ~SimplePolynomial() {};

    T& operator[](int index) const {
        return this->koeffs[index];
    }

    virtual const T operator()(const T point) const {
        T mult = 1;
        return koeffs.template reduce<T>(0, [&](const T acc, int i, const Point<T>& vector) -> const T {
            T result = mult * vector[i];
            mult *=point;
            return acc + result;
        });
    }

    virtual const T operator()(const Point<T>& point) const {
        return this->operator()(point[0]);
    }

    /* change this.koeffs **/
    virtual const SimplePolynomial<T> operator()(const SimplePolynomial<T>& polynomial) const {
        SimplePolynomial<T> accPolynomial(1, 1);
        return SimplePolynomial<T>(koeffs.template reduce<Point<T>>(Point<T>(getDegree() * polynomial.getDegree() + 1, 0.0), [&](const Point<T> acc, int i, const Point<T>& a) -> Point<T> {
            accPolynomial.koeffs.forEach([&](int j, const Point<T>& b) -> void {
                acc[j]+= b[j] * a[i];
            });
            accPolynomial *=polynomial;
            return acc;
        }));
    }

    SimplePolynomial<T>& replaceVariable(const SimplePolynomial<T>& polynomial) {
        SimplePolynomial<T> accPolynomial(1, 1);
        koeffs = koeffs.template reduce<Point<T>>(Point<T>(getDegree() * polynomial.getDegree() + 1, 0.0), [&](const Point<T> acc, int i, const Point<T>& a) -> Point<T> {
            accPolynomial.koeffs.forEach([&](int j, const Point<T>& b) -> void {
                acc[j]+= b[j] * a[i];
            });
            accPolynomial *=polynomial;
            return acc;
        });
        return *this;
    }

    SimplePolynomial<T>& derivative() {
        if(koeffs.length() <=1) {
            koeffs = Point<T>(1, 0.0);
        }
        Point<T> res(koeffs.length() - 1);
        res.forEach([&](int i, const Point<T>& point) -> void {
           res[i] = koeffs[i+1]*(i+1);
        });
        koeffs = res;
        return *this;
    }

    const SimplePolynomial<T> getDerivative() const {
        if(koeffs.length() <=1) {
            return SimplePolynomial<T>(1, 0.0);
        }
        Point<T> res(koeffs.length() - 1);
        res.forEach([&](int i, const Point<T>& point) -> void {
            res[i] = koeffs[i+1]*(i+1);
        });
        return SimplePolynomial<T>(res);
    }

    SimplePolynomial<T>& operator=(const SimplePolynomial<T>& polynomial) {
        koeffs = polynomial.getKoeffs();
        return *this;
    }

    const SimplePolynomial<T> operator+(const SimplePolynomial<T>& polynomial) const {
        if(getDegree() >= polynomial.getDegree()){
            return SimplePolynomial<T>(this->getKoeffs() + polynomial.getKoeffs());
        }
        return SimplePolynomial<T>(polynomial.getKoeffs() + this->getKoeffs());
    }

    SimplePolynomial<T>& operator+=(const SimplePolynomial<T>& polynomial) {
        if(getDegree() >= polynomial.getDegree()) {
            koeffs += polynomial.getKoeffs();
        } else {
            koeffs = polynomial.getKoeffs() + koeffs;
        }
        return *this;
    }

    const SimplePolynomial<T> operator-(const SimplePolynomial<T>& polynomial) const {
        if(getDegree() >= polynomial.getDegree()){
            return SimplePolynomial<T>(this->getKoeffs() - polynomial.getKoeffs());
        }
        return SimplePolynomial<T>(polynomial.getKoeffs() - this->getKoeffs());
    }

    SimplePolynomial<T>& operator-=(const SimplePolynomial<T>& polynomial) {
        if(getDegree() >= polynomial.getDegree()) {
            koeffs -= polynomial.getKoeffs();
        } else {
            koeffs = polynomial.getKoeffs() - koeffs;
        }
        return *this;
    }

    const SimplePolynomial<T> operator*(const T scalar) const {
        return SimplePolynomial<T>(this->getKoeffs() * scalar);
    }

    SimplePolynomial<T>& operator*=(const T scalar) {
        koeffs *= scalar;
        return *this;
    }

    const SimplePolynomial<T> operator*(const SimplePolynomial<T>& polynomial) const {
        return SimplePolynomial<T>(koeffs.template reduce<Point<T>>(Point<T>(getDegree() + polynomial.getDegree() + 1, 0.0), [&](Point<T> acc, int i, const Point<T>& a) -> Point<T> {
            polynomial.koeffs.forEach([&](int j, const Point<T>& b) -> void {
                acc[i+j]+= a[i] * b[j];
            });
            return acc;
        }));
    }
    SimplePolynomial<T>& operator*=(const SimplePolynomial<T>& polynomial) {
        koeffs = koeffs.template reduce<Point<T>>(Point<T>(getDegree() + polynomial.getDegree() + 1, 0.0), [&](Point<T> acc, int i, const Point<T>& a) -> Point<T> {
            polynomial.koeffs.forEach([&](int j, const Point<T>& b) -> void {
                acc[i+j]+= a[i] * b[j];
            });
            return acc;
        });
        return *this;
    }



    void print(std::string name, int countOfNumbers=0, int precision=0) const {
        std::cout<<std::endl<<name<<std::endl;
        koeffs.forEachReverse([&](int i, const Point<T>& point) -> void {
            if(i != point.length() - 1) {
                std::cout<<" "<<(sign(point[i]) ? "+" : "");
            }
            if(countOfNumbers != 0 && precision != 0) {
                std::cout<<std::fixed<<std::setw(countOfNumbers)<<std::setprecision(precision)<<point[i];
            } else {
                std::cout<<point[i];
            }
            if(i > 1) {
                std::cout<<"x^("<<i<<")";
            } else if(i == 1) {
                std::cout<<"x";
            }
        });
        std::cout<<std::endl;
    }

    const Point<T> getKoeffs() const {
        return Point<T>(this->koeffs);
    }

    const int getDegree() const {
        return koeffs.length() - 1;
    }

};

void TestSimplePolynomial();


#endif //SIMPLE_POLYNOMIAL_H