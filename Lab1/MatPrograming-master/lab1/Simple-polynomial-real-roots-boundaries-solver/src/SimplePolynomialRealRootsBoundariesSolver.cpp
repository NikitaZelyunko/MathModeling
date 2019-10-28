

#include <SimplePolynomialRealRootsBoundariesSolver.h>

template <class T>
void printBoundaries(Point<Point<T>> boundaries) {
    for(int i = 0; i < 2; i++) {
        boundaries[i].print(i ? "negative boundaries:" : "positive boundaries:");
    }
}

void TestSimplePolynomialRealRootsBoundariesSolver() {
    std::cout<<std::endl<<std::endl;
    std::cout<<"--------------------------"<<std::endl;
    std::cout<<"Test: Simple polynomial real roots boundaries solver test."<<std::endl;
    SimplePolynomial<double> polynomial1(3,4);
    polynomial1.print("f1(x):");
    SimplePolynomialRealRootsBoundariesSolver<double> solver1(polynomial1, 0, 0.00001);
    printBoundaries(solver1.solve());

    Point<double> koeffs2(4);
    koeffs2[0] = 8; koeffs2[1] = 2; koeffs2[2]= -5; koeffs2[3] = 1;
    SimplePolynomial<double> polynomial2(koeffs2);
    polynomial2.print("f2(x):");
    SimplePolynomialRealRootsBoundariesSolver<double> solver2(polynomial2, 0, 0.00001);
    printBoundaries(solver2.solve());

    polynomial2.print("f2(x):");
    SimplePolynomialRealRootsBoundariesSolver<double> solver3(polynomial2, 1, 0.00001);
    printBoundaries(solver3.solve());

    Point<double> koeffs4(7);
    koeffs4[0] = -2880; koeffs4[1] = -1344; koeffs4[2] = 1780; koeffs4[3] = 72; koeffs4[4] = -161; koeffs4[5] = 12; koeffs4[6] = 1;
    SimplePolynomial<double> polynomial4(koeffs4);
    polynomial4.print("f4(x):");
    SimplePolynomialRealRootsBoundariesSolver<double> solver4(koeffs4.length(), Point<double>::toArray(koeffs4), 1, 1);
    printBoundaries(solver4.solve());

//    polynomial4.print("f4(x):");
//    SimplePolynomialRealRootsBoundariesSolver<double> solver5(polynomial4, 1, 0.00001);
//    printBoundaries(solver5.solve());

    std::cout<<"--------------------------"<<std::endl;
}

