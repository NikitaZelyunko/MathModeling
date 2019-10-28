#include "../include/Utils.h"

double eNorm(double * x, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++) {
		result += pow(x[i], 2);
	}
	return sqrt(result);
}

void printArray(double * a, int n, string comment)
{
	if (a) {
		cout << endl << comment << endl;
		for (int i = 0; i < n; i++) {
			cout << a[i] << " ";
		}
		cout << endl;
	}
	else cout << "printArray: a not exist" << endl;
}

void printArray(double ** a, int n, int m, string comment)
{
	if (a) {
		cout << comment << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				cout << a[i][j] << " ";
			}
			cout << endl;
		}
	}
	else cout << "printArray: a not exist" << endl;
}

double ** createTwoDimArray(int n, int m)
{
	double** result = new double*[n];
	for (int i = 0; i < n; i++) {
		result[i] = new double[n];
	}
	return result;
}

void deleteTwoDimArray(double ** a, int n)
{
	for (int i = 0; i < n; i++) {
		delete[] a[i];
	}
	delete a;
}

double * getCopyArray(double * a, int n)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		result[i] = a[i];
	}
	return result;
}

double ** getCopyArray(double ** a, int n, int m)
{
	double **result = createTwoDimArray(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result[i][j] = a[i][j];
		}
	}
	return result;
}

double * multOnScalar(double * x, int n, double scalar)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		result[i] = scalar*x[i];
	}
	return result;
}

double ** multOnScalar(double ** A, int n, int m, double scalar)
{
	double** result = createTwoDimArray(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result[i][j] = scalar * A[i][j];
		}
	}
	return result;
}

double * expVector(double * x, double exp, int n)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		result[i] = pow(x[i], exp);
	}
	return result;
}

double * sumVectors(double * x1, double* x2, int n)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		result[i] = x1[i] + x2[i];
	}
	return result;
}

double * sumVectorsMultOnScalars(double a, double * x1, double b, double * x2, int n)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		result[i] = a*x1[i] + b*x2[i];
	}
	return result;
}

double * substrVectors(double * x1, double * x2, int n)
{
	double* minusX2 = multOnScalar(x2, n, -1);
	for (int i = 0; i < n; i++) {
		minusX2[i] += x1[i];
	}
	return minusX2;
}

double scalarMult(double * x1, double* x2, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++) {
		result += x1[i] * x2[i];
	}
	return result;
}

double * multMatrixOnVector(double ** A, int n, int m, double * x)
{
	double* result = createOneDimArray(n);
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < m; j++) {
			sum += A[i][j] * x[j];
		}
		result[i] = sum;
	}
	return result;
}

double ** MakeIdentityMatrix(int n)
{
	double** E = createTwoDimArray(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			E[i][j] = 0;
		}
		E[i][i] = 1;
	}
	return E;
}