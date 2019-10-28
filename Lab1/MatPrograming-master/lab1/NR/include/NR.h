#ifndef NR_H
#define NR_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <limits>

#include "Utils.h"

using namespace std;

//swap k and l rows
void swapRows(double** A, int m, int k, int l);
int findRowIndexMaxElementInColumnFromL(double ** A, int n, int k, int l);


class MyCustomOneDimFunc {
	F func;
	Grad grad;
	Hessian hessian;

private: 

	double* xK = nullptr;
	double* dK = nullptr;
	int dimCount;

	double** dKMatrix = nullptr;

	inline double getPointsInInterval(double a, double b, double h) {
		return floor((b - a) / h) + 1;
	}

	inline double* getArg(double t) {
		return sumVectorsMultOnScalars(1, xK, t, dK, dimCount);
	}

	double** getDkMatrix() {
		double** result = createTwoDimArray(dimCount, dimCount);
		for (int i = 0; i < dimCount; i++) {
			for (int j = 0; j < dimCount; j++) {
				result[i][j] = dK[i] * dK[j];
			}
		}
		return result;
	}

public:
	MyCustomOneDimFunc(F func, Grad grad, Hessian hessian, double* xK, double* dK, int dimCount) {
		this->func = func;
		this->grad = grad;
		this->hessian = hessian;
		this->xK = xK;
		this->dK = dK;
		this->dimCount = dimCount;
	}

	~MyCustomOneDimFunc() {
		// if (dKMatrix) {
		// 	deleteTwoDimArray(dKMatrix, dimCount);
		// }
	}

	double operator()(double t) {
		double* arg = getArg(t);
		double result = func(arg);
		delete arg;
		return result;
	}

	double proizv1(double t) {
		double* arg = getArg(t);
		double* gradient = grad(arg, dimCount);
		double result = scalarMult(gradient, dK, dimCount);
		delete arg;
		delete gradient;
		return result;
	}

	// double proizv2(double t) {
	// 	double* arg = getArg(t);
	// 	double** hessian = this->hessian(arg, dimCount);
	// 	if (!dKMatrix) {
	// 		dKMatrix = getDkMatrix();
	// 	}
	// 	double result = 0;
	// 	for (int i = 0; i < dimCount; i++) {
	// 		for (int j = 0; j < dimCount; j++) {
	// 			result += hessian[i][j] * dKMatrix[i][j];
	// 		}
	// 	}
	// 	deleteTwoDimArray(hessian, dimCount);
	// }

	double getGlobalMinimum(double a, double b, double epsilon) {
		double h = epsilon;
		// count of point in interval
		double n = getPointsInInterval(a, b, h);
		cout << "n:" << n << endl;
		vector<vector<double>> extremumsIntervals;
		double y0 = proizv1(a);
		double y1 = 0;
		bool sY0 = sign(y0);
		bool sY1;
		double x0 = a;
		double x1 = a + h;
		for (int i = 1; i < n; i++, x0 = x1, x1 +=h, y0 = y1, sY0 = sY1) {
			y1 = proizv1(x1);
			sY1 = sign(y1);
			if (sY0 != sY1) {
				cout << "interval:" << "[" << x0 << "," << x1 << "]" << endl;
				vector<double> extrInterval;
				extrInterval.push_back(x0);
				extrInterval.push_back(x1);
				extremumsIntervals.push_back(extrInterval);
			}
		}
		return 0;
	}
};

class MatrixWithDeterm {
public:
	double** matrix;
	int n;
	int m;
	double determ;

	MatrixWithDeterm(double** matrix, int n, int m, double determ) {
		this->matrix = matrix;
		this->n = n;
		this->m = m;
		this->determ = determ;
	}

	void print(string comment) {
		if (matrix) {
			cout << comment << endl;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					cout << fixed << setw(7)<< setprecision(3)<<  matrix[i][j] << " ";
				}
				cout << endl;
			}
		}
		else cout << "printArray: " << comment << " not exist"<<endl;
	}

	~MatrixWithDeterm() {
		if (matrix != nullptr) {
			for (int i = 0; i < n; i++)
				delete[] matrix[i];
			delete matrix;
		}
	}
};

MatrixWithDeterm* inverseMatrix(MatrixWithDeterm &matrix, double epsilon);

#endif