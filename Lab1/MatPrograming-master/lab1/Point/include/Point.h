#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <string>
#include <iomanip>

#include <functional>

template<class T>
class Point {
    private:
        int n = 0;
        T* coeffs;
    // dimension


    bool checkBoundary(int index) const {
        if(index >= 0) {
            if(index < n) {
                return true;
            }
        }
        return false;
    }

    inline int correctDimension(int dimension) const {
        return dimension < 0 ? 0 : dimension;
    }

    inline void printOutOfRange(int index) const {
        std::cout<<"Out of range "<<"n="<<n<<" index="<<index<<std::endl;
    }

    inline void fillCoeffs(const T& filler) {
        for(int i = 0; i < n; i++)
            this->coeffs[i] = filler;
    }

    inline T* createOneDimArray(int n) const {
        return new T[n];
    }

    public:

    static T* copyArray(T* x, int n) {
        T* result = new T[n];
        for(int i = 0; i < n; i++) 
            result[i] = x[i];
        return result;
    }

    static T* toArray(const Point<T>& point){
//        T* arr = new T[point.length()];
        return point.reduce<T*>(new T[point.length()], [](T* acc, int i, const Point<T> point) -> T* {
           acc[i] = T(point[i]);
           return acc;
        });
    }

    Point() {
        this->n = 0;
        this->coeffs = nullptr;
    }

    Point(int n) {
        this->n = correctDimension(n);
        this->coeffs = createOneDimArray(this->n);
    }

    Point(int n, const T& filler) {
        this->n = correctDimension(n);
        this->coeffs = new T[this->n];
        fillCoeffs(filler);
    }

    Point(int n, T* coeffs) {
        this->n = correctDimension(n);
        this->coeffs = copyArray(coeffs, this->n);
    }

    Point(const Point<T>& x) {
        this->n = x.n;
        this->coeffs = copyArray(x.coeffs, x.n);
    }

    ~Point() {
        delete coeffs;
    }

    int length() const {
        return n;
    }

    T& operator[] (int index) const {
        if(checkBoundary(index)) {
            return coeffs[index];
        }
        printOutOfRange(index);
        return coeffs[n-1];
    }

    inline T& enorm() const {
        T result = 0;
        for(int i =0; i < n; i++)
            result+=coeffs[i];
        return sqrt(result);
    }

    inline Point<T>& operator =(const Point<T>& x) {
        n = x.n;
        delete coeffs;
        coeffs = copyArray(x.coeffs, x.n);
        return *this;
    }

    inline const Point<T> operator =(const T& scalar) const {
        delete coeffs;
        fillCoeffs(scalar);
        return *this;
    }

    inline const Point<T> operator *(const T& scalar) const {
        Point<T> res(*this);
        for(int i = 0; i < n; i++)
            res[i] *=scalar;
        return res;
    }

    inline Point<T>& operator *=(const T& scalar) {
        for(int i = 0; i < n; i++)
            coeffs[i]*=scalar;
        return *this;
    }    
    
    inline const Point<T> operator /(const T& scalar) const{
        Point<T> res(*this);
        if(scalar != 0) {
            for(int i = 0; i < n; i++)
                res[i] /=scalar;
        }
        return res;
    }

    inline const Point<T> operator /=(const T& scalar) {
        if(scalar != 0) {
            for(int i = 0; i < n; i++)
                coeffs[i]/=scalar;
        }
        return Point(*this);
    }

    inline const Point<T> operator +(const Point<T>& x) const {
        Point<T> res(*this);
        for(int i = 0; i < x.length(); i++)
            res[i] += x[i];
        return res;
    }

    inline Point<T>& operator +=(const Point<T>& x) {
        for(int i = 0; i < x.length(); i++)
            coeffs[i]+=x[i];
        return *this;
    }

    inline const Point<T> operator -(const Point<T>& x) const {
        Point<T> res(*this);
        for(int i = 0; i < x.length(); i++)
            res[i] -= x[i];
        return res;
    }

    inline const Point<T> operator -=(const Point<T>& x) {
        for(int i = 0; i < x.length(); i++)
            coeffs[i]-=x[i];
        return Point(*this);
    }

    inline void print(std::string name, int countOfNumbers=0, int precision=0) const {
        std::cout<<std::endl<<name<<std::endl;
        this->forEach([&](int i, const Point<T>& point) -> void {
            if(countOfNumbers != 0 && precision != 0) {
                std::cout<<std::fixed<<std::setw(countOfNumbers)<<std::setprecision(precision)<<point[i]<<" ";
            } else {
                std::cout<<point[i]<<" ";
            }

    });
        std::cout<<std::endl;
    }

    inline void forEach(std::function<void (int, const Point<T>&)> callback) const  {
        for(int i = 0; i < n; i++) {
            callback(i, *this);
        }
    }

    inline void forEachReverse(std::function<void (int, const Point<T>&)> callback) const  {
        for(int i = n-1; i >= 0; i--) {
            callback(i, *this);
        }
    }

    inline void forEachFrom(int from, std::function<void (int, const Point<T>&)> callback) const  {
        for(int i = from; i < n; i++) {
            callback(i, *this);
        }
    }

    inline void forEachFromTo(int from, int to, std::function<void (int, const Point<T>&)> callback) const  {
        for(int i = from; i < to && i < n; i++) {
            callback(i, *this);
        }
    }

    // todo try to add this to matrix
    template<class R>
    inline const R reduce(const R startAcc, std::function<const R (const R&, int, const Point<T>&)> callback) const {
        R bufAcc = startAcc;
        for(int i = 0; i < n; i++){
            bufAcc = callback(bufAcc, i, *this);
        }
        return bufAcc;
    }

    template<class R>
    inline const R reduceReverse(const R startAcc, std::function<const R (const R&, int, const Point<T>&)> callback) const {
        R bufAcc = startAcc;
        for(int i = n-1; i >=0; i--){
            bufAcc = callback(bufAcc, i, *this);
        }
        return bufAcc;
    }

};

void TestPoint();

#endif
