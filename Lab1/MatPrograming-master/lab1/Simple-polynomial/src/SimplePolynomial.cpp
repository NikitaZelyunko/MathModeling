#include "SimplePolynomial.h"

void TestSimplePolynomial() {
    std::cout<<std::endl<<std::endl;
    std::cout<<"--------------------------"<<std::endl;
    std::cout<<"Test: Simple polynomial test."<<std::endl;
    SimplePolynomial<double> polinom(3, 2);
    polinom.print("f(x):");
    std::cout<<"f(2):"<<polinom(2)<<std::endl;

    Point<double> koeffs1(4, 3);
    SimplePolynomial<double> polinom1(koeffs1);
    polinom1.print("f1(x):");
    std::cout<<"f1(3):"<<polinom1(3)<<std::endl;

    double* koeffs2 = new double[4];
    koeffs2[0] = 2; koeffs2[1] = 1.5; koeffs2[2] = -4; koeffs2[3] = 0;
    SimplePolynomial<double> polinom2(4, koeffs2);
    polinom2.print("f2(x):");
    std::cout<<"f2(1.5):"<<polinom2(1.5)<<std::endl;

    SimplePolynomial<double> polinom3(polinom2);
    polinom3.print("f3(x):");
    std::cout<<"f3(1.5):"<<polinom3(1.5)<<std::endl;

    Point<double> point1(2, 2);
    std::cout<<"f3(2,2):"<<polinom3(point1)<<std::endl;

    polinom3.getKoeffs().print("f3.koeffs:");

    (polinom + polinom2).print("f(x) + f2(x):");
    (polinom2 + polinom).print("f2(x) + f(x):");
    (polinom + polinom).print("f(x) + f(x):");

    polinom += polinom2;
    polinom.print("f(x) + f2(x):");

    (polinom1 - polinom2).print("f1(x) - f2(x):");
    polinom1 = polinom2;
    polinom1.print("f1(x) = f2(x):");
    polinom1 = polinom1 * 2;
    polinom1.print("f1(x) = 2 * f1(x):");
    polinom1 *= 2;
    polinom1.print("f1(x) *= 2:");

    Point<double> koeffs4(3);
    koeffs4[0] = -3; koeffs4[1] = 4; koeffs4[2] = 2;
    SimplePolynomial<double> polynomial4(koeffs4);
    polynomial4.print("f4(x):");
    Point<double> koeffs5(4);
    koeffs5[0] = 2; koeffs5[1] = -3; koeffs5[2] = 4; koeffs5[3] = 1;
    SimplePolynomial<double> polynomial5(koeffs5);
    polynomial5.print("f5(x):");
    (polynomial4*polynomial5).print("f4(x) * f5(x):");
    (polynomial5*polynomial4).print("f5(x) * f4(x):");
    polynomial4 *=polynomial5;
    polynomial4.print("f4(x) *= f5(x)");
    Point<double> koeffs6(2);
    koeffs6[0] = 3; koeffs6[1] = 2;
    SimplePolynomial<double> polynomial6(koeffs6);
    polynomial6.print("f6(x):");
    polynomial5(polynomial6).print("f5(f6(x))");
    polynomial6.replaceVariable(polynomial5).print("f6(f5(x))");
    polynomial6[0] = 11;
    polynomial6.print("f6(x)");
    polynomial6.derivative().print("f6'(x):");
    std::cout<<"--------------------------"<<std::endl;
}
