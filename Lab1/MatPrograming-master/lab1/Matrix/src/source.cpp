#include "../include/Matrix.h"

void TestMatrix() {
    cout<<endl<<endl;
    cout<<"-------------------"<<endl;
    cout<<"----Test Matrix----"<<endl;
    Matrix<double> a = Matrix<double>(2,2,-1.5);
    a.print("a:");

    a[0][0] = 0;
    a[1][1] = 10;
    a.print("a:");

    double** coeffs = new double*[3];
    for(int i = 0; i < 3; i++) {
        coeffs[i] = new double[2];
        for(int j = 0; j < 2; j++)
            coeffs[i][j] = 1;
    }

    Matrix<double> b = Matrix<double>(3, 2, coeffs);
    b.print("b:");

    Matrix<double> identity = Matrix<double>::createIdentityMatrix(4);
    identity.print("E:");

    Matrix<double> c = Matrix<double>(identity);

    c[0][0] = -1;
    c[1][1] = -1;
    c[2][2] = -1;
    c[3][3] = -1;
    c.print("c:");

    Matrix<double> d = Matrix<double>(2,2,2.0);
    Matrix<double> e = Matrix<double>::createIdentityMatrix(2);
    d.print("d before:");
    e.print("e before:");
    d-=e;
    d.print("d after:");
    e.print("e after:");

    Matrix<double> f = d*e;
    f.print("d*e:");

    e *=d;

    e.print("e after after:");

    Matrix<double> f1 = Matrix<double>(2, 3);
    f1[0][0] = 1; f1[0][1] = 2; f1[0][2] = 3;
    f1[1][0] = 4; f1[1][1] = 5; f1[1][2] = 6;

    Matrix<double> f2 = Matrix<double>(3,2);
    f2[0][0] = 1; f2[0][1] = 2;
    f2[1][0] = 3; f2[1][1] = 4;
    f2[2][0] = 5; f2[2][1] = 6;

    f1.print("f1:");
    f2.print("f2:");

    Matrix<double> f3 = f1*f2;
    f3.print("f1 * f2:");

    f1*=f2;
    f1.print("f1*f2:");

    Matrix<double> o1 = Matrix<double>::createIdentityMatrix(3);
    o1.print("o1:");
    o1+=o1;
    o1.print("o1+o1:");
    o1-= o1 * 2;
    o1.print("o1+01 - 2*(o1+o1):");
    o1*=3;
    o1.print("-//-*3:");
    o1 = o1*3 - o1;
    o1.print("o1*3 - o1:");
    Matrix<double> o2 = Matrix<double>(3,3,4);
    o2.print("o2:");
    o1 = o2 + o1;
    o1.print("o2 + o1:");

    f2.transpose().print("f2.transpose():");

    Matrix<double>::transpose(f2).print("f2.transpose again:");
    f2.print("f2:");
    cout<<"--End test Matrix--"<<endl;
    cout<<"-------------------"<<endl;
    cout<<endl<<endl;
}