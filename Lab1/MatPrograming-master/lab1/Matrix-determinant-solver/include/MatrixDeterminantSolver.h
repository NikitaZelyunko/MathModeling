#ifndef MATRIX_DETERMINANT_SOLVER
#define MATRIX_DETERMINANT_SOLVER

#include <iostream>
#include <string>

#include "MatrixAbstractSolver.h"
#include "Matrix.h"

const int MATRIX_IS_DEGENERATE = 0;

template <class R>
class MatrixDeterminantSolver : public MatrixAbstractSolver<R, R>
{
private:

    const R computeTriangularDet(const Matrix<R>& x) const {
        if(Matrix<R>::isSquare(x)) {
            return x.forDiagReduce(1, [&](R acc, int i, const Matrix<R>& matrix) -> const R {
                return acc * matrix[i][i];
            });
        }
        this->matrixNotSquare();
        return this->getResult();
    }

    const R evaluateGeneralFormMatrixDet() const {
        R detMult = 1;
        Matrix<R> resultMatrix = getLowerTriangleMatrix(this->matrix, detMult, this->getEpsilon());
        if(Matrix<R>::isSquare(resultMatrix)) {
            return detMult * computeTriangularDet(resultMatrix);
        }
        else {
            this->matrixNotSquare();
        }
        return this->getResult();
    }

protected: 
    virtual const R computeResult() const {
        if(this->isTriangular) {
            return this->computeTriangularDet(this->matrix);
        } 
        else {
            return this->evaluateGeneralFormMatrixDet();
        }
    }

public:
    using MatrixAbstractSolver<R, R>::MatrixAbstractSolver;

    void print(string prefix) const {
        if(this->matrix) {
            cout<<"Determinant matrix:"<<endl;
            this->matrix.print(prefix);
            cout<<"Is: "<<this->getResult()<<endl;
        } else {
            this->solverNotInited();
        }
    }

    static int findRowIndexMaxElementInColumnFromL(const Matrix<R>& x, int k, int l) {
        int indexMax = l;
        R maxElem = x[l][l];
        for (int i = l; i < x.getRowCount(); i++) {
            if (abs(x[i][k]) > abs(maxElem)) {
                indexMax = i;
                maxElem = x[i][k];
            }
        }
        return indexMax;
    }

    static const Matrix<R> getLowerTriangleMatrix(const Matrix<R>& matrix, R& detMult, R epsilon) {
        Matrix<R> resultMatrix = Matrix<R>(matrix);
        if(Matrix<R>::isSquare(resultMatrix)) {
            try {
                resultMatrix.forDiag([&](int i, const Matrix<R>& x) -> void {
                    if(i != x.getRowCount() - 1) {
                        int maxIndex = findRowIndexMaxElementInColumnFromL(x, i, i);
                        if (i != maxIndex) {
                            x.swapRows(i, maxIndex);
                            detMult *= -1;
                        }
                        R xii = x[i][i];
                        if (abs(xii) < epsilon) {
                            throw MATRIX_IS_DEGENERATE;
                        }
                        for(int k = i + 1; k < x.getColumnCount(); k++) {
                            R identifier = x[k][i]/xii;
                            for (int j = i; j < x.getColumnCount(); j++) {
                                x[k][j] -= x[i][j] * identifier;
                            }
                        }
                    }
                });
            }
            catch (int err) {
                if(err == MATRIX_IS_DEGENERATE) {
                    detMult = 0;
                }
            }
        }
        return resultMatrix;
    }
};

void TestMatrixDeterminantSolver();

#endif