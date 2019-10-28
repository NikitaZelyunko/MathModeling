#ifndef MATRIX_SYLVESTER_TEST_H
#define MATRIX_SYLVESTER_TEST_H

#include "Point.h"
#include "Matrix.h"
#include "MatrixAbstractSolver.h"

const int ALTERNATING_CERTAINTY = 0;
const int POSITIVE_CERTAINTY = 1;
const int NEGATIVE_CERTAINTY = 2;

template <class S>
class MatrixSylvesterTest : public MatrixAbstractSolver<int, S> {
    private:
    const Point<S> resolveCertainty() const {
        int minDimension = min(this->matrix.getColumnCount(), this->matrix.getRowCount());
        Point<S> cornerMinors = Point<S>(minDimension, 1);
        Matrix<S> resultMatrix = Matrix<S>(this->matrix);
        if(Matrix<S>::isSquare(resultMatrix)) {
            try {
                resultMatrix.forDiag([&](int i, const Matrix<S>& x) -> void {
                    S xii = x[i][i];
                    cornerMinors[i] *= xii;
                    if(abs(xii) >= this->getEpsilon()) {
                        if(i + 1 < cornerMinors.length()) {
                            cornerMinors[i+1] *= cornerMinors[i];
                        }
                    } else {
                        int notNullRow = this->findFirstNotNullElementIndexInColumn(x, i, i);
                        if(notNullRow != i) {
                            cornerMinors.forEachFromTo(i + 1, notNullRow, [&](int i, const Point<S> &point) -> void {
                                point[i] = 0;
                            });
                            x.swapRows(i, notNullRow);
                            xii = x[i][i];
                            cornerMinors[notNullRow] *= -xii;

                        } else {
                            cornerMinors.forEachFrom(i+1, [&](int i, const Point<S> &point) -> void {
                                point[i] = 0;
                            });
                            throw ALTERNATING_CERTAINTY;
                        }
                    }
                    for(int k = i + 1; k < x.getColumnCount(); k++) {
                        S identifier = x[k][i]/xii;
                        for (int j = i; j < x.getColumnCount(); j++) {
                            x[k][j] -= x[i][j] * identifier;
                        }
                    }
                });
            }
            catch (int err) {}
        } else {
            return Point<S>(minDimension, 0.0);
        }
        return cornerMinors;
    }

    protected:

    virtual const int computeResult() const {
        Point<S> cornerMinors = resolveCertainty();
        cornerMinors.print("cornerMinors:");
        try {
            return cornerMinors.template reduce<int>(0, [&](const int& acc, int i, const Point<S>& point) -> const int& {
                if(i == 0) {
                    if(point[i] >= this->getEpsilon()) {
                        return POSITIVE_CERTAINTY;
                    }
                    if(point[i] <= -this->getEpsilon()) {
                        return NEGATIVE_CERTAINTY;
                    }
                }
                if(acc == POSITIVE_CERTAINTY && point[i] >= this->getEpsilon()) {
                    return acc;
                }
                if(acc == NEGATIVE_CERTAINTY) {
                    if(i % 2 == 0) {
                        if(point[i] <= -this->getEpsilon()) {
                            return acc;
                        }
                    } else {
                        if(point[i] >= this->getEpsilon()) {
                            return acc;
                        }
                    }
                }
                throw ALTERNATING_CERTAINTY;
            });
        }
        catch (int err) {
            return ALTERNATING_CERTAINTY;
        }
    }

    public:
    using MatrixAbstractSolver<int, S>::MatrixAbstractSolver;
    /*
        If matrix is positive certainty return 1;
        If matrix is negative certainty return 2;
        Else return 0.
    **/
    virtual int solve() {
        return MatrixAbstractSolver<int, S>::solve();
    }
};

void TestMatrixSylvesterTest();

#endif
