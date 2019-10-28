#include "../include/Point.h"
#include "Utils.h"

void TestPoint() {
    std::cout<<std::endl<<std::endl;
    std::cout<<"--------------"<<std::endl;
    std::cout<<"TEST POINT LIB"<<std::endl;
    Point<double> a = Point<double>(1);
    Point<double> b = Point<double>(2);
    Point<double> c = Point<double>(2, -2);
    Point<double> d = Point<double>(2, 0.0);
    double *coeffs1 = createOneDimArray(3);
    for(int i = 0; i < 3; i++) {
        coeffs1[i] = i;
    }
    Point<double> e = Point<double>(3, coeffs1);

    Point<double> f = Point<double>(e);

    a.print("a:");
    b.print("b:");
    c.print("c:");
    d.print("d:");
    e.print("e:");
    f.print("f:");

    Point<double> g = f*3;
    g.print("g:");

    a[0] = -1;

    b[0] = 6;
    b[1] = -7;
    a.print("a:");
    b.print("b:");

    a = b;
    a.print("a:");
    cout<<"--------------"<<endl;
}
