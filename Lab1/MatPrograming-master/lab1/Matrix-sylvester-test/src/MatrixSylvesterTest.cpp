#include "../include/MatrixSylvesterTest.h"


void TestMatrixSylvesterTest() {
    cout<<endl<<endl;
    cout<<"--------------------"<<endl;
    cout<<"Test: Matrix Sylvester test."<<endl;
    Matrix<double> a(3,3);
    a[0][0] = 1; a[0][1] = 2; a[0][2] = 3;
    a[1][0] = 4; a[1][1] = 2; a[1][2] = 2;
    a[2][0] = 2; a[2][1] = 9; a[2][2] = 7;
    a.print("a:");
    MatrixSylvesterTest<double> tester1(a, 0.00001);
    cout<<"Result:"<<tester1.solve()<<endl;

    Matrix<double> b(3,3);
    b[0][0] = 0; b[0][1] = 2; b[0][2] = 0;
    b[1][0] = 2; b[1][1] = 0; b[1][2] = 0;
    b[2][0] = 0; b[2][1] = 0; b[2][2] = 0;
    b.print("b:");
    MatrixSylvesterTest<double> tester2(b, 0.00001);
    cout<<"Result:"<<tester2.solve()<<endl;

    Matrix<double> c(3,3);
    c[0][0] = -2; c[0][1] = 2; c[0][2] = 0;
    c[1][0] = 2; c[1][1] = 3; c[1][2] = 0;
    c[2][0] = 0; c[2][1] = 0; c[2][2] = 2;
    c.print("c:");
    MatrixSylvesterTest<double> tester3(c, 0.00001);
    cout<<"Result:"<<tester3.solve()<<endl;

    Matrix<double> d(3,3);
    d[0][0] = -2; d[0][1] = 2; d[0][2] = 0;
    d[1][0] = 2; d[1][1] = -3; d[1][2] = 0;
    d[2][0] = 0; d[2][1] = 0; d[2][2] = -2;
    d.print("d:");
    MatrixSylvesterTest<double> tester4(d, 0.00001);

    Matrix<double> e(2,2);
    e[0][0] = 2; e[0][1] = 0;
    e[1][0] = 0; e[1][1] = 2;
    e.print("e:");
    MatrixSylvesterTest<double> tester5(e, 0.00001);
    cout<<"Result:"<<tester5.solve()<<endl;

    Matrix<double> f(4,4);
    f[0][0] = 2; f[0][1] = 0; f[0][2] = 1; f[0][3] = 3;
    f[1][0] = 7; f[1][1] = 2; f[1][2] = 4; f[1][3] = 1;
    f[2][0] = 2; f[2][1] = 2; f[2][2] = 2; f[2][3] = 0;
    f[3][0] = -1; f[3][1] = 2; f[3][2] = 3; f[3][3] = 2;
    f.print("f:");
    MatrixSylvesterTest<double> tester6(f, 0.00001);
    cout<<"Result:"<<tester6.solve()<<endl;

    Matrix<double> g(4,4,0.0);
    g.print("g:");
    MatrixSylvesterTest<double> tester7(g, 0.00001);
    cout<<"Result:"<<tester7.solve()<<endl;

    Matrix<double> h = Matrix<double>::createIdentityMatrix(4);
    h.print("h:");
    MatrixSylvesterTest<double> tester8(h, 0.00001);
    cout<<"Result:"<<tester8.solve()<<endl;

    Matrix<double> i = Matrix<double>::createIdentityMatrix(4) * (-1);
    i.print("i:");
    MatrixSylvesterTest<double> tester9(i, 0.00001);
    cout<<"Result:"<<tester9.solve()<<endl;
    cout<<"--------------------"<<endl;
}