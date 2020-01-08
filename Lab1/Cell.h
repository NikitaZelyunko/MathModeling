//
// Created by nikit on 14.10.2019.
//

#ifndef LAB1_CELL_H
#define LAB1_CELL_H

#include "iostream"
#include "vector"

#include "Point.h"

const int IN = 0;
const int OUT = 1;
const int WALL = 2;
const int INNER = 3;

template <class T>
class Cell : public Point<T> {
protected:
    int type;

    bool checkType(int type) const {
        return type >= 0 && type <= 3;
    }

    std::vector<Cell<T>*> neighbors;
public:

    using Point<T>::Point;

    Cell(): Point<T>() {
        this->type = INNER;
    }

    Cell(int n): Point<T>(n) {
        this->type = INNER;
    }

    Cell(int n, const T& filler): Point<T>(n, filler) {
        this->type = INNER;
    }

    Cell(int n, const T& filler, int type): Cell<T>(n, filler) {
        if(checkType(type)) {
            this->type = type;
        }
    }

    Cell(int n, T* coeffs): Point<T>(n, coeffs) {
        this->type = INNER;
    }

    Cell(int n, T* coeffs, int type): Cell(n, coeffs) {
        if(checkType(type)) {
            this->type = type;
        }
    }

    Cell(const Cell<T>& x): Point<T>(x) {
        this->type = x.type;
    }

    void setType(int type) {
        this->type = type;
    }

    int getType() {
        return this->type;
    }

    bool isInner() const {
        return  type == INNER;
    }

    bool isIn() const {
        return type == IN;
    }

    bool  isOut() const {
        return type == OUT;
    }

    bool  isWall() const {
        return  type == WALL;
    }

    bool isBorder() const {
        return isIn() || isOut() || isWall();
    }

    inline Cell<T>& operator =(const Point<T>& x) {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator=(x);
        }
        return (*this);
    }

    inline const Cell<T> operator =(const T& scalar) const {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator=(scalar);
        }
        return (*this);
    }

    inline Cell<T>& operator *=(const T& scalar) {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator*=(scalar);
        }
        return (*this);
    }

    inline const Cell<T> operator /=(const T& scalar) {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator/=(scalar);
        }
        return (*this);
    }

    inline Cell<T>& operator +=(const Point<T>& x) {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator+=(x);
        }
        return (*this);
    }

    inline const Cell<T> operator -=(const Point<T>& x) {
        if(this->isInner() || this->isOut() || this->isWall()) {
            Point<T>::operator-=(x);
        }
        return (*this);
    }

    void addNeighbor(Cell<T>* neighbor) {
        this->neighbors.push_back(neighbor);
    }

    Cell<T>& getNeighbor(int index) {
        return *(this->neighbors[index]);
    }
};

#endif //LAB1_CELL_H
