#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "functional"

#include "Cell.h"

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

template <class T>
void print2DArray(T** arr, int n, int m) {
    for(int i=0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            std::cout<<arr[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

template <class T>
T*** copyArray(T*** arr, long int n, long int m, long int l) {
    T*** res = create3DArray<T>(n,m,l);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < l; k++) {
                res[i][j][k] = arr[i][j][k];
            }
        }
    }
    return res;
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

	T h = (x1-x0) / N;

	T x = x0;
	T y = y0;
	T z = z0;
	
	T t = t0;

	long int nx = N;
	long int ny = N;
	long int nz = N;

	T*** curLayer = create3DArray<T>(nx, ny, nz);
	T*** nextLayer = create3DArray<T>(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                curLayer[i][j][k] = 0;
            }
        }
    }

	for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k] = 1;
            nextLayer[0][j][k] = 1;
            curLayer[nx-1][j][k] = 1;
            nextLayer[nx-1][j][k] = 1;
        }
    }

    T h_2 = h*h;
    T a = tau / h_2;
	for (int time = 0; time < M; time++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
				    T res = curLayer[i + 1][j][k] + curLayer[i - 1][j][k] +
                            curLayer[i][j + 1][k] + curLayer[i][j - 1][k] +
                            curLayer[i][j][k + 1] + curLayer[i][j][k - 1] -
                            6 * curLayer[i][j][k];
				    T ijk = curLayer[i][j][k];
					nextLayer[i][j][k] = a*(
						curLayer[i + 1][j][k] + curLayer[i - 1][j][k] +
						curLayer[i][j + 1][k] + curLayer[i][j - 1][k] +
						curLayer[i][j][k + 1] + curLayer[i][j][k - 1] -
						6 * curLayer[i][j][k]) + curLayer[i][j][k];
					std::cout<<"";
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

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                curLayer[i][j][k] = 0;
            }
        }
    }

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k] = 1;
            bufLayer[0][j][k] = 1;
            curLayer[nx-1][j][k] = 1;
            bufLayer[nx-1][j][k] = 1;
        }
    }

    T h_2 = h*h;
    T koeff = tau/h_2;

    for (int time = 0; time < M; time++) {

        // F
        T F = 0;

        for (int i = 0; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    F = (curLayer[i+1][j][k] - curLayer[i][j][k]) * koeff;
                    if(i == 0) {
                        bufLayer[i+1][j][k]-=F;
                    } else if(i == nx-2) {
                        bufLayer[i][j][k]+=F;
                    } else {
                        bufLayer[i+1][j][k]-=F;
                        bufLayer[i][j][k]+=F;
                    }
                }
            }
        }

        // G

        for (int j = 0; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                for (int k = 1; k < nz - 1; k++) {
                    F = (curLayer[i][j+1][k] - curLayer[i][j][k]) * koeff;
                    if(j == 0) {
                        bufLayer[i][j+1][k]-=F;
                    } else if(j == ny - 2) {
                        bufLayer[i][j][k]+=F;
                    } else {
                        bufLayer[i][j+1][k]-=F;
                        bufLayer[i][j][k]+=F;
                    }
                }
            }
        }

        // H

        for (int k = 0; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    F = (curLayer[i][j][k+1] - curLayer[i][j][k]) * koeff;
                    if(k == 0) {
                        bufLayer[i][j][k+1]-=F;
                    } else if(k == nz - 2) {
                        bufLayer[i][j][k]+=F;
                    } else {
                        bufLayer[i][j][k+1]-=F;
                        bufLayer[i][j][k]+=F;
                    }
                }
            }
        }
        delete3DArray(curLayer, nx, ny);
        curLayer = copyArray(bufLayer, nx, ny, nz);
    }

    delete3DArray<T>(bufLayer, nx, ny);
    return curLayer;
}

template<class T>
void abstractIterateMethod(
        Cell<T>*** cells,
        long int nx, long int ny, long int nz,
        std::function<void(const Cell<T>&, long int, long int, long int, Cell<T>***)> callback
        ) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                callback(cells[i][j][k], i, j, k, cells);
            }
        }
    }
}
template<class T>
T* arrSum(const T& a, T* arr1, const T& a1, T* arr2, const T& a2, long int n) {
    T* res = create1DArray<T>(n);
    for(int i = 0; i < n; i++) {
        res[i] = a*(a1 * arr1[i] + a2*arr2[i]);
    }
    return res;
}

template<class T>
T* accumulateSum(T* acc, const T& a, T* arr, int n) {
    for(int i = 0; i < n; i++){
        acc[i]+=a*arr[i];
    }
    return acc;
}

template<class T>
void gasDynamic(Cell<T>*** cells,
                long int nx, long int ny, long int nz,
                T hx, T hy, T hz,
                T tau, T t0, T t1,
                long int paramCount
        ) {
    T hx_2 = hx * hx;
    T koeffX = tau/hx_2;
    T hy_2 = hy * hy;
    T koeffY = tau/hy_2;
    T hz_2 = hz * hz;
    T koeffZ = tau/hz_2;

    long int M = floor((t1 - t0) / tau + 1);
    std::cout<<"All bad"<<std::endl;
    for (int time = 0; time < M; time++) {
        std::cout<<"All bad:"<<time<<std::endl;
        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++) {
                for(int k = 0; k < nz; k++) {
                    Cell<T>& cell = cells[i][j][k];
//                    T* nullParams = create1DArray<T>(paramCount);
//                    Point<Cell<T>*>xNeighbors = Point<Cell<T>*>(2);
//                    Point<Cell<T>*>yNeighbors = Point<Cell<T>*>(2);
//                    Point<Cell<T>*>zNeighbors = Point<Cell<T>*>(2);
//                    if(!cell.isBorder()) {
//                        xNeighbors[0] = &cells[i-1][j][k];
//                        xNeighbors[1] = &cells[i+1][j][k];
//
//                        yNeighbors[0] = &cells[i][j-1][k];
//                        yNeighbors[1] = &cells[i][j+1][k];
//
//                        zNeighbors[0] = &cells[i][j][k-1];
//                        zNeighbors[1] = &cells[i][j][k+1];
//                    } else {
//                        if(i == 0) {
//                            xNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                            xNeighbors[1] = &cells[i+1][j][k];
//                        } else if(i == nx-1) {
//                            xNeighbors[0] = &cells[i-1][j][k];
//                            xNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                        } else {
//                            xNeighbors[0] = &cells[i-1][j][k];
//                            xNeighbors[1] = &cells[i+1][j][k];
//                        }
//
//                        if(j == 0) {
//                            yNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                            yNeighbors[1] = &cells[i][j+1][k];
//                        } else if(j == ny - 1) {
//                            yNeighbors[0] = &cells[i][j-1][k];
//                            yNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                        } else {
//                            yNeighbors[0] = &cells[i][j-1][k];
//                            yNeighbors[1] = &cells[i][j+1][k];
//                        }
//
//                        if(k == 0) {
//                            zNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                            zNeighbors[1] = &cells[i][j][k+1];
//                        } else if(k == nz - 1) {
//                            zNeighbors[0] = &cells[i][j][k-1];
//                            zNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                        } else {
//                            zNeighbors[0] = &cells[i][j][k-1];
//                            zNeighbors[1] = &cells[i][j][k+1];
//                        }
//                    }

                    Point<T> F1 = (cell - cell.getNeighbor(0))*koeffX;
                    Point<T> F2 = (cell.getNeighbor(1) - cell)*koeffX;

                    Point<T> G1 = (cell - cell.getNeighbor(2))*koeffY;
                    Point<T> G2 = (cell.getNeighbor(3) - cell)*koeffY;

                    Point<T> H1 = (cell - cell.getNeighbor(4))*koeffZ;
                    Point<T> H2 = (cell.getNeighbor(5) - cell)*koeffZ;

                    cell.getNeighbor(0)+=F1;
                    cell.getNeighbor(1)-=F2;

                    cell.getNeighbor(2)+=G1;
                    cell.getNeighbor(3)-=G2;

                    cell.getNeighbor(4)+=H1;
                    cell.getNeighbor(5)-=H2;
                }
            }
        }
//        abstractIterateMethod<T>(
//                cells,
//                nx, ny, nz,
//                [&](const Cell<T>& cell, long int i, long int j, long int k, Cell<T>*** cells) -> void {
//                T* nullParams = create1DArray<T>(paramCount);
//                Point<Cell<T>*>xNeighbors = Point<Cell<T>*>(2);
//                Point<Cell<T>*>yNeighbors = Point<Cell<T>*>(2);
//                Point<Cell<T>*>zNeighbors = Point<Cell<T>*>(2);
//                if(!cell.isBorder()) {
//                    xNeighbors[0] = &cells[i-1][j][k];
//                    xNeighbors[1] = &cells[i+1][j][k];
//
//                    yNeighbors[0] = &cells[i][j-1][k];
//                    yNeighbors[1] = &cells[i][j+1][k];
//
//                    zNeighbors[0] = &cells[i][j][k-1];
//                    zNeighbors[1] = &cells[i][j][k+1];
//                } else {
//                    if(i == 0) {
//                        xNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                        xNeighbors[1] = &cells[i+1][j][k];
//                    } else if(i == nx-1) {
//                        xNeighbors[0] = &cells[i-1][j][k];
//                        xNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                    } else {
//                        xNeighbors[0] = &cells[i-1][j][k];
//                        xNeighbors[1] = &cells[i+1][j][k];
//                    }
//
//                    if(j == 0) {
//                        yNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                        yNeighbors[1] = &cells[i][j+1][k];
//                    } else if(j == ny - 1) {
//                        yNeighbors[0] = &cells[i][j-1][k];
//                        yNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                    } else {
//                        yNeighbors[0] = &cells[i][j-1][k];
//                        yNeighbors[1] = &cells[i][j+1][k];
//                    }
//
//                    if(k == 0) {
//                        zNeighbors[0] = new Cell<T>(paramCount, nullParams);
//                        zNeighbors[1] = &cells[i][j][k+1];
//                    } else if(k == nz - 1) {
//                        zNeighbors[0] = &cells[i][j][k-1];
//                        zNeighbors[1] = new Cell<T>(paramCount, nullParams);
//                    } else {
//                        zNeighbors[0] = &cells[i][j][k-1];
//                        zNeighbors[1] = &cells[i][j][k+1];
//                    }
//                }
//
//                Point<T> F1 = (cell - *(xNeighbors[0]))*koeffX;
//                Point<T> F2 = (*(xNeighbors[1]) - cell)*koeffX;
//
//                Point<T> G1 = (cell - *(yNeighbors[0]))*koeffY;
//                Point<T> G2 = (*(yNeighbors[1]) - cell)*koeffY;
//
//                Point<T> H1 = (cell - *(zNeighbors[0]))*koeffZ;
//                Point<T> H2 = (*(zNeighbors[1]) - cell)*koeffZ;
//
//                *(xNeighbors[0])+=F1;
//                *(xNeighbors[1])-=F2;
//
//                *(yNeighbors[0])+=G1;
//                *(yNeighbors[1])-=G2;
//
//                *(zNeighbors[0])+=H1;
//                *(zNeighbors[1])-=H2;
//                delete1DArray(nullParams);
//        });
    }
}
