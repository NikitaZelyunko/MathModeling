#include "../include/MatrixDeterminantSolver.h"

void TestMatrixDeterminantSolver() {
    Matrix<double> a(2,2,2);
    a.print("a:");
    MatrixDeterminantSolver<double> solver = MatrixDeterminantSolver<double>(a, 0.00001, true);
    cout<<"isSolved:"<<solver.isSolved()<<endl;
    cout<<"Solve result:"<<solver.solve()<<endl;
    cout<<"isSolved:"<<solver.isSolved()<<endl;

    MatrixDeterminantSolver<double> solver2 = MatrixDeterminantSolver<double>(a, 0.00001);
    cout<<"isSolved:"<<solver2.isSolved()<<endl;
    cout<<"Solve result:"<<solver2.solve()<<endl;
    cout<<"isSolved:"<<solver2.isSolved()<<endl;

    Matrix<double> i = Matrix<double>::createIdentityMatrix(2);
    MatrixDeterminantSolver<double> solver3 = MatrixDeterminantSolver<double>(i, 0.00001);
    cout<<"isSolved:"<<solver3.isSolved()<<endl;
    cout<<"Solve result:"<<solver3.solve()<<endl;
    cout<<"isSolved:"<<solver3.isSolved()<<endl;

    Matrix<double> b(3,3,0.0);
    b[0][0] = 1; b[0][1] = 2;  b[0][2] = 3;
    b[1][0] = 4; b[1][1] = 2;  b[1][2] = 2;
    b[2][0] = 2; b[2][1] = 9;  b[2][2] = 7;
    b.print("b:");
    MatrixDeterminantSolver<double> solver4 = MatrixDeterminantSolver<double>(b, 0.00001);
    cout<<"isSolved:"<<solver4.isSolved()<<endl;
    cout<<"Solve result:"<<solver4.solve()<<endl;
    cout<<"isSolved:"<<solver4.isSolved()<<endl;

    Matrix<double> c(1,1,222);
    c.print("c:");
    MatrixDeterminantSolver<double> solver5 = MatrixDeterminantSolver<double>(c, 0.00001);
    cout<<"isSolved:"<<solver5.isSolved()<<endl;
    cout<<"Solve result:"<<solver5.solve()<<endl;
    cout<<"isSolved:"<<solver5.isSolved()<<endl;
}