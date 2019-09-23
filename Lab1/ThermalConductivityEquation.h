#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

template<class T>
T* create1DArray(int n) {
	return new T[n];
}

template<class T>
T** create2DArray(int n, int m) {
	T** arr = create1DArray<T*>(n);
	for (int i = 0; i < n; i++)
		arr[i] = create1DArray<T>(m);
	return arr;
}

template<class T>
T*** create3DArray(int n, int m, int l) {
	T*** arr = create2DArray<T*>(n, m);
	for (int i = 0; i < n; i++) 
		for(int j = 0; j < m; j++)
			arr[i][j] = create1DArray<T>(m);
	return arr;
}

template<class T>
void delete1DArray(T* arr) {
	delete arr;
}

template<class T>
void delete2DArray(T** arr, int n) {
	for (int i = 0; i < n; i++)
		delete1DArray(arr[i]);
	delete arr;
}

template<class T>
void delete3DArray(T*** arr, int n, int m) {
	for (int i = 0; i < n; i++)
		delete2DArray(arr[i], m);
	delete arr;
}

template<class T>
T*** Reshenie_Uravn_Teploprovodnosti_methodom_progonki_yavn(
	long int N,
	T x0, T x1,
	T y0, T y1,
	T z0, T z1,

	T tau, T t0, T t1
)
{
	long int M = floor(t1 - t0) / tau + 1;

	T h = 1.0 / N;

	T x = x0;
	T y = y0;
	T z = z0;
	
	T t = t0;

	long int nx = N;
	long int ny = N;
	long int nz = N;

	T*** curLayer = create3DArray<T>(nx, ny, nz);
	T*** nextLayer = create3DArray<T>(nx, ny, nz);

	for (int j = 0; j < ny; j++) {
		for (int k = 0; k < nz; k++) {
			curLayer[0][j][k] = 1;
			nextLayer[0][j][k] = 1;
		}
	}

	for (int i = 1; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				curLayer[i][j][k] = 0;
				nextLayer[i][j][k] = 0;
			}
		}
	}

	T h_2 = h*h;
	T a = tau * h_2;

	for (int time = 0; time < M; time++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
					nextLayer[i][j][k] = a*(
						curLayer[i + 1][j][k] + curLayer[i - 1][j][k] +
						curLayer[i][j + 1][k] + curLayer[i][j - 1][k] +
						curLayer[i][j][k + 1] + curLayer[i][j][k - 1] -
						6 * curLayer[i][j][k]) + curLayer[i][j][k];
				}
			}
		}
		T*** buf = curLayer;
		curLayer = nextLayer;
		nextLayer = buf;
	}

	delete3DArray<T>(nextLayer, nx, ny);
	return curLayer;
}

template<class T>
T*** Reshenie_Uravn_Teploprovodnosti_flux(
        long int N,
        T x0, T x1,
        T y0, T y1,
        T z0, T z1,

        T tau, T t0, T t1
)
{
    long int M = floor(t1 - t0) / tau + 1;

    T h = 1.0 / N;

    T x = x0;
    T y = y0;
    T z = z0;

    T t = t0;

    long int nx = N;
    long int ny = N;
    long int nz = N;

    T*** curLayer = create3DArray<T>(nx, ny, nz);
    T*** bufLayer = create3DArray<T>(nx, ny, nz);

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k] = 1;
            bufLayer[0][j][k] = 1;
        }
    }

    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                curLayer[i][j][k] = 0;
                bufLayer[i][j][k] = 0;
            }
        }
    }

    T h_2 = h*h;
    T a = tau * h_2;

    for (int time = 0; time < M; time++) {
        std::cout<<t+tau*time<<std::endl;

        // F
        T F = 0;

        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                F = (curLayer[1][j][k] - curLayer[0][j][k])/h;
                bufLayer[1][j][k]-=F;
            }
        }

        for (int i = 1; i < nx - 2; i++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    F = (curLayer[i+1][j][k] - curLayer[i][j][k])/h;
                    bufLayer[i+1][j][k]-=F;
                    bufLayer[i][j][k]+=F;
                }
            }
        }

        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                F = (curLayer[nx-1][j][k] - curLayer[nx-2][j][k])/h;
                bufLayer[nx-2][j][k]+=F;
            }
        }

        // G

        for (int i = 1; i < nx - 1; i++) {
            for (int k = 1; k < nz - 1; k++) {
                F = (curLayer[i][1][k] - curLayer[i][0][k])/h;
                bufLayer[i][1][k]-=F;
            }
        }

        for (int j = 1; j < ny - 2; j++) {
            for (int i = 1; i < nx - 1; i++) {
                for (int k = 1; k < nz - 1; k++) {
                    F = (curLayer[i][j+1][k] - curLayer[i][j][k])/h;
                    bufLayer[i][j+1][k]-=F;
                    bufLayer[i][j][k]+=F;
                }
            }
        }
        for (int i = 1; i < nx - 1; i++) {
            for (int k = 1; k < nz - 1; k++) {
                F = (curLayer[i][ny-1][k] - curLayer[i][ny-2][k])/h;
                bufLayer[i][ny-2][k]+=F;
            }
        }

        // H
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                F = (curLayer[i][j][1] - curLayer[i][j][0])/h;
                bufLayer[i][j][1]-=F;
            }
        }

        for (int k = 1; k < nz - 2; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    F = (curLayer[i][j][k+1] - curLayer[i][j][k])/h;
                    bufLayer[i][j][k+1]-=F;
                    bufLayer[i][j][k]+=F;
                }
            }
        }

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                F = (curLayer[i][j][nz-1] - curLayer[i][j][nz-2])/h;
                bufLayer[i][j][nz-2]+=F;
            }
        }

        T koeff = tau/h;
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    bufLayer[i][j][k]*=koeff;
                }
            }
        }
        T*** swapBuf = curLayer;
        curLayer = bufLayer;
        bufLayer = swapBuf;
    }

    delete3DArray<T>(bufLayer, nx, ny);
    return curLayer;
}