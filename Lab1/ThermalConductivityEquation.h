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
void thermalDynamic(Cell<T>*** cells,
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
    }
}

template<class T>
const Point<T> getF(const Cell<Point<T>>& cell) {
    Point<T> F = Point<T>(cell[0].length());
    F[0] = cell[0][1];
    F[1] = pow(F[0],2)/cell[0][0] + cell[0][5];
    F[2] = F[0]*cell[0][2]/cell[0][0];
    F[3] = F[0]*cell[0][3]/cell[0][0];
    F[4] = (cell[0][4] + cell[0][5])*F[0]/cell[0][0];
    F[5] = 0;
    return F;
}

template<class T>
const Point<T> getG(const Cell<Point<T>>& cell) {
    Point<T> G = Point<T>(cell[0].length());
    G[0] = cell[0][2];
    G[1] = G[0]*cell[0][1]/cell[0][0];
    G[2] = pow(G[0], 2)/cell[0][0] + cell[0][5];
    G[3] = G[0]*cell[0][3]/cell[0][0];
    G[4] = (cell[0][4] + cell[0][5])*G[0]/cell[0][0];
    G[5] = 0;
    return G;
}

template<class T>
const Point<T> getH(const Cell<Point<T>>& cell) {
    Point<T> H = Point<T>(cell[0].length());
    H[0] = cell[0][3];
    H[1] = H[0] * cell[0][1]/cell[0][0];
    H[2] = H[0] * cell[0][2] / cell[0][0];
    H[3] = pow(H[0], 2)/cell[0][0] + cell[0][5];
    H[4] = (cell[0][4] + cell[0][5])*H[0]/cell[0][0];
    H[5] = 0;
    return H;
}

template<class T>
T getAlpha(Point<T> firstCell, Point<T> secondCell) {
    T firstC = sqrt(firstCell[5]*1.4/firstCell[0]);
    T secondC = sqrt(secondCell[5]*1.4/secondCell[0]);

    T firstU = firstCell[1]/firstCell[0];
    T firstV = firstCell[2]/firstCell[0];
    T firstW = firstCell[3]/firstCell[0];

    T secondU = secondCell[1]/secondCell[0];
    T secondV = secondCell[2]/secondCell[0];
    T secondW = secondCell[3]/secondCell[0];

    return std::max(sqrt(pow(firstU, 2) + pow(firstV, 2) + pow(firstW, 2)) + firstC, sqrt(pow(secondU, 2) + pow(secondV, 2) + pow(secondW, 2)) + secondC);
}


template<class T>
void gasDynamic(Cell<Point<T>>*** cells,
                long int nx, long int ny, long int nz,
                T hx, T hy, T hz,
                T tau, T t0, T t1,
                int paramCount
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
                    Cell<Point<T>>& cell = cells[i][j][k];

                    // F
                    Point<T> Fcur = getF(cell);
                    Point<T> Fleft = getF(cell.getNeighbor(0));
                    Point<T> Fright = getF(cell.getNeighbor(1));

                    T FalphaLeft = getAlpha(cell.getNeighbor(0)[0], cell[0]);
                    T FalphaRight = getAlpha(cell[0], cell.getNeighbor(1)[0]);
                    // G
                    Point<T> Gcur = getG(cell);
                    Point<T> Gleft = getG(cell.getNeighbor(2));
                    Point<T> Gright = getG(cell.getNeighbor(3));

                    T GalphaLeft = getAlpha(cell.getNeighbor(2)[0], cell[0]);
                    T GalphaRight = getAlpha(cell[0], cell.getNeighbor(3)[0]);
                    // H
                    Point<T> Hcur = getH(cell);
                    Point<T> Hleft = getH(cell.getNeighbor(4));
                    Point<T> Hright = getH(cell.getNeighbor(5));

                    T HalphaLeft = getAlpha(cell.getNeighbor(4)[0], cell[0]);
                    T HalphaRight = getAlpha(cell[0], cell.getNeighbor(5)[0]);

                    Point<T> F_prev_2 = (Fcur + Fleft - ((cell[0] - cell.getNeighbor(0)[0]) * FalphaLeft))*0.5;
                    Point<T> F_next_2 = (Fright + Fcur - (cell.getNeighbor(1)[0] - cell[0]) * FalphaRight)*0.5;

                    Point<T> G_prev_2 = (Gcur +  Gleft - (cell[0] - cell.getNeighbor(2)[0]) * GalphaLeft)*0.5;
                    Point<T> G_next_2 = (Gright + Gcur - (cell.getNeighbor(3)[0] - cell[0]) * GalphaRight)*0.5;

                    Point<T> H_prev_2 = (Hcur + Hleft - (cell[0] - cell.getNeighbor(4)[0])* HalphaLeft)*0.5;
                    Point<T> H_next_2 = (Hright + Hcur - (cell.getNeighbor(5)[0] - cell[0]) * HalphaRight)*0.5;

                    // U

//                    Point<T> F1 = (cell - cell.getNeighbor(0))*koeffX;
//                    Point<T> F2 = (cell.getNeighbor(1) - cell)*koeffX;
//
//                    Point<T> G1 = (cell - cell.getNeighbor(2))*koeffY;
//                    Point<T> G2 = (cell.getNeighbor(3) - cell)*koeffY;
//
//                    Point<T> H1 = (cell - cell.getNeighbor(4))*koeffZ;
//                    Point<T> H2 = (cell.getNeighbor(5) - cell)*koeffZ;
//
//                    cell.getNeighbor(0)+=F1;
//                    cell.getNeighbor(1)-=F2;
//
//                    cell.getNeighbor(2)+=G1;
//                    cell.getNeighbor(3)-=G2;
//
//                    cell.getNeighbor(4)+=H1;
//                    cell.getNeighbor(5)-=H2;
                }
            }
        }
    }
}
