#ifndef MATRIX_ABSTRACT_SOLVER
#define MATRIX_ABSTRACT_SOLVER

#include "AbstractSolver.h"
#include "Matrix.h"

template<class R, class S>
class MatrixAbstractSolver : public AbstractSolver<R, S> {
    protected:
        Matrix<S> matrix = Matrix<S>(0,0);
        bool isTriangular = false;
    
        virtual void matrixNotSquare() const {
            std::cout<<std::endl<<"The matrix is not square. The search for a determinant cannot be performed."<<std::endl;
        }

        bool checkTriangular(const Matrix<S>& matrix, bool cache) {
            bool isTriang = isAnyTriangular(matrix);
            if(cache) {
                this->isTriangular = isTriang;
            }
            return isTriang;
        }

    public:
        using AbstractSolver<R,S>::AbstractSolver;

        MatrixAbstractSolver(const Matrix<S>& matrix, S epsilon): AbstractSolver<R, S>(epsilon) {
            this->matrix = Matrix<S>(matrix);
        }

        MatrixAbstractSolver(const Matrix<S>& matrix, S epsilon, bool isTriangular): MatrixAbstractSolver<R, S>(matrix, epsilon) {
            this->isTriangular = isTriangular;
        }

        bool isLowerTriangular(const Matrix<S>& matrix) const {
            return matrix.everyLowerTriangle([&](int i, int j, const Matrix<S>& matrix) -> bool {
                return abs(matrix[i][j]) <= this->getEpsilon();
            });
        }

        bool isUpperTriangular(const Matrix<S>& matrix) const {
            return matrix.everyUpperTriangle([&](int i, int j, const Matrix<S>& matrix) -> bool {
                return abs(matrix[i][j]) <= this->getEpsilon();
            });
        }

        bool isAnyTriangular(const Matrix<S>& matrix) const {
            return isUpperTriangular(matrix) || isLowerTriangular(matrix);
        }

        int findFirstNotNullElementIndexInColumn(const Matrix<S>& matrix, int columnIndex, int startFind) const {
            for(int i = startFind + 1; i < matrix.getRowCount(); i++) {
                if(abs(matrix[i][columnIndex]) >= this->getEpsilon()) {
                    return i;
                }
            }
            return startFind;
        }
};

#endif