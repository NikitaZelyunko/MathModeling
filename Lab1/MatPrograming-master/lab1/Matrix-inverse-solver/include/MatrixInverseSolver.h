#ifndef MATRIX_INVERSE_SOLVER
#define MATRIX_INVERSE_SOLVER

#include "Matrix.h"
#include "MatrixAbstractSolver.h"
#include "MatrixDeterminantSolver.h"

template<class S>
class MatrixInverseSolver : public MatrixAbstractSolver<Matrix<S>, S> {

    protected:
    virtual const Matrix<S> computeResult() const {
        S detMult = 1;
        Matrix<S> resultMatrix = Matrix<S>(this->matrix);
        Matrix<S> E = Matrix<S>::createIdentityMatrix(resultMatrix.getRowCount());

        if(Matrix<S>::isSquare(resultMatrix)) {
            try {
                resultMatrix.forDiag([&](int i, const Matrix<S>& x) -> void {
                    int maxIndex = MatrixDeterminantSolver<S>::findRowIndexMaxElementInColumnFromL(x, i, i);
                    if (i != maxIndex) {
                        x.swapRows(i, maxIndex);
                        E.swapRows(i, maxIndex);
                        detMult *= -1;
                    }
                    S xii = x[i][i];
                    if (abs(xii) < this->getEpsilon()) {
                        throw MATRIX_IS_DEGENERATE;
                    }
                    for(int j = i; j < x.getColumnCount(); j++) {
                        x[i][j] /= xii;
                    }
                    for (int j = 0; j < x.getColumnCount(); j++) {
                        E[i][j] /= xii;
                    }
                    for(int k = i + 1; k < x.getColumnCount(); k++) {
                        S identifier = x[k][i];
                        for (int j = i; j < x.getColumnCount(); j++) {
                            x[k][j] -= x[i][j] * identifier;
                            E[k][j] -= E[i][j] * identifier;
                        }
                    }
                });
            }
            catch (int err) {
                if(err == MATRIX_IS_DEGENERATE) {
                    detMult = 0;
                }
            }

            MatrixDeterminantSolver<S> detSolver = MatrixDeterminantSolver<S>(resultMatrix, this->getEpsilon(), true);
            S det = detMult * detSolver.solve();
            if(abs(det) >= this->getEpsilon()) {
                for (int i = resultMatrix.getRowCount() - 1; i >= 0; i--) {
                    for (int k = i - 1; k >= 0; k--) {
                        S aki = resultMatrix[k][i];
                        // A[k][i] = 0;
                        for (int j = 0; j < resultMatrix.getColumnCount(); j++) {
                            E[k][j] -= E[i][j] * aki;
                        }
                    }
                }
            }
        }
        return E;
    }

    public:
        using MatrixAbstractSolver<Matrix<S>, S>::MatrixAbstractSolver;
};

void TestMatrixInverseSolver();
#endif