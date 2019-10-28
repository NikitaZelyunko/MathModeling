#include "../include/Header.h"

double f3(double* x) {
	return pow(pow(x[0],2) + x[1] -11, 2) + pow(x[0] + pow(x[1], 2) - 7, 2);
}

double* gradF3(double* x, int n) {
	double* result = new double[n];

	result[0] = 4 * pow(x[0], 3) + 4 * x[0] * x[1] - 42 * x[0] + 2 * pow(x[1], 2) - 14;
	result[1] = 4 * pow(x[1], 3) + 2 * pow(x[0], 2) + 4 * x[0] * x[1] - 26 * x[1] - 22;
	return result;
}

double ** gessianF3(double* x, int n) {
	double** result = new double*[n];
	for(int i = 0; i < n; i++) {
		result[i] = new double[n];
	}
	result[0][0] = 12 * pow(x[0], 2) + 4 * x[1] - 42; result[0][1] = 4 * (x[0] + x[1]);
	result[1][0] = 4 * (x[0] + x[1]); result[1][1] = 12 * pow(x[1], 2) + 4 * x[0] - 26;
	
	return result;
}

double f1(double* x) {
	return pow(x[0] - 2, 2) + pow(x[1] - 5, 2);
}

double* gradF1(double* x, int n) {
	double* result = new double[n];
	result[0] = 2 * (x[0] - 2);
	result[1] = 2 * (x[1] - 5);
	return result;
}

double** gessianF1(double* x, int n) {
	double** result = new double*[n];
	for(int i = 0; i < n; i++) {
		result[i] = new double[n];
	}
	result[0][0] = 2; result[0][1] = 0;
	result[1][0] = 0; result[1][1] = 2;
	return result;
}

double fXtD1(double* x, double* d, double t) {
	return pow(x[0] + t*d[0] - 2, 2) + pow(x[1] + t * d[1] - 5,2);
}

double dFxtD1(double* x, double* d, double t) {
	return 2 * (d[0] * (x[0] + t* d[0] - 2) + d[1] * (x[1] + t * d[1] - 5));
}

double d2FxtD1(double* x, double* d, double t) {
	return 2 * (pow(d[0], 2) + pow(d[1], 2));
}

double getMinimum1(double* x, double* d) {
	return - (d[0] * (x[0] - 2) + d[1] * (x[1] - 5)) / (pow(d[0], 2) + pow(d[1], 2));
}


double f7(double* x) {
	return pow(x[0], 3) + pow(x[1], 2) - x[0] * (x[1] - 2) + 3 * x[1] - 4;
}

double* gradF7(double* x, int n) {
	double* result = new double[n];
	result[0] = 3 * pow(x[0], 2) - x[1] - 2;
	result[1] = - x[0] + 2 * x[1] + 3;
	return result;
}

double** gessianF7(double* x, int n) {
	double** result = new double*[n];
	for (int i = 0; i < n; i++) {
		result[i] = new double[n];
	}
	result[0][0] = 6 * x[0]; result[0][1] = -1;
	result[1][0] = -1; result[1][1] = 2;
	return result;
}


int main() {
	int n = 2;

	double *x0 = new double[n];
	x0[0] = 0;
	x0[1] = 0;
//	TestMatrix();
//	TestMatrixDeterminantSolver();
//	TestMatrixInverseSolver();
//    TestSimplePolynomial();
    TestSimplePolynomialRealRootsBoundariesSolver();
//    TestMatrixSylvesterTest();
//    double* result = getTochnValue(x0, n, 0.001, 0.001, f1, gradF1, gessianF1, getMinimum1, 1000000);
//    MyCustomOneDimFunc myFunc(f1, gradF1, gessianF1, x0, x0, n);
//    printArray(result, n, "Result:");
	return 0;
}