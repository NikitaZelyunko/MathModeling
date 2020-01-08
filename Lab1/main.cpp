
#include "header.h"

void printCSV(double *** arr,
              long int N,
              double x0, double x1,
              double y0, double y1,
              double z0, double z1,
              std::string filename
) {

    double h = 1.0 / N;

    double x = x0;
    double y = y0;
    double  z = z0;

    long int nx = N;
    long int ny = N;
    long int nz = N;

    std::ofstream file(filename.c_str());
    file << "x,y,z,T" << std::endl;
    for (int i = 0; i < nx; i++, x+=h) {
        y = y0;
        z = z0;
        for (int j = 0; j < ny; j++, y+=h) {
            z = z0;
            for (int k = 0; k < nz; k++, z+=h) {
                file<< x << "," << y << "," << z << "," <<arr[i][j][k]<< std::endl;
            }
        }
    }
    file.close();
}

template <class T>
void solveNeighbors(
        Cell<T>*** curLayer,
        Cell<T>& nullCell,
        long int nx, long int ny, long int nz,
        long int paramCount
) {
    for(int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Cell<T>& cell = curLayer[i][j][k];
                if (!cell.isBorder()) {
                    cell.addNeighbor(&(curLayer[i - 1][j][k]));
                    cell.addNeighbor(&(curLayer[i + 1][j][k]));

                    cell.addNeighbor(&(curLayer[i][j - 1][k]));
                    cell.addNeighbor(&(curLayer[i][j + 1][k]));

                    cell.addNeighbor(&(curLayer[i][j][k - 1]));
                    cell.addNeighbor(&(curLayer[i][j][k + 1]));
                } else {
                    if (i == 0) {
                        cell.addNeighbor(&nullCell);
                        cell.addNeighbor(&(curLayer[i + 1][j][k]));
                    } else if (i == nx - 1) {
                        cell.addNeighbor(&(curLayer[i - 1][j][k]));
                        cell.addNeighbor(&nullCell);
                    } else {
                        cell.addNeighbor(&(curLayer[i - 1][j][k]));
                        cell.addNeighbor(&(curLayer[i + 1][j][k]));
                    }

                    if (j == 0) {
                        cell.addNeighbor(&nullCell);
                        cell.addNeighbor(&(curLayer[i][j + 1][k]));
                    } else if (j == ny - 1) {
                        cell.addNeighbor(&(curLayer[i][j - 1][k]));
                        cell.addNeighbor(&nullCell);
                    } else {
                        cell.addNeighbor(&(curLayer[i][j - 1][k]));
                        cell.addNeighbor(&(curLayer[i][j + 1][k]));
                    }

                    if (k == 0) {
                        cell.addNeighbor(&nullCell);
                        cell.addNeighbor(&(curLayer[i][j][k + 1]));
                    } else if (k == nz - 1) {
                        cell.addNeighbor(&(curLayer[i][j][k - 1]));
                        cell.addNeighbor(&nullCell);
                    } else {
                        Cell<T>& zero = curLayer[0][0][0];
                        cell.addNeighbor(&(curLayer[i][j][k - 1]));
                        cell.addNeighbor(&(curLayer[i][j][k + 1]));
                    }
                }
            }
        }
    }
}

template <class T>
  Cell<T>*** generateThermalCells(
        long int nx, long int ny, long int nz,
        T hx, T hy, T hz,
        T tau, T t0, T t1,
        long int paramCount
        ) {


    Cell<T>*** curLayer = create3DArray<Cell<T>>(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                curLayer[i][j][k] = Cell<T>(paramCount, 0.0, 3);
            }
        }
    }

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k] = Cell<T>(paramCount, 1, 0);
            curLayer[nx-1][j][k] = Cell<T>(paramCount, 1, 0);
        }
    }

    for(int i = 1; i < nx-1; i++) {
        for(int k = 0; k < nz; k++){
            curLayer[i][0][k] = Cell<T>(paramCount, 0.0, 0);
            curLayer[i][ny-1][k] = Cell<T>(paramCount, 0.0, 0);
        }
    }
    for(int i = 1; i < nx - 1; i++) {
        for(int j = 0; j < ny; j++){
            curLayer[i][j][0] = Cell<T>(paramCount, 0.0, 0);
            curLayer[i][j][nz-1] = Cell<T>(paramCount, 0.0, 0);
        }
    }

    Cell<T>& nullCell = *(new Cell<T>(paramCount, 0.0, 1));
      for(int i = 0; i < nx; i++) {
          for (int j = 0; j < ny; j++) {
              for (int k = 0; k < nz; k++) {
                  Cell<T>& cell = curLayer[i][j][k];
                  if (!cell.isBorder()) {
                      cell.addNeighbor(&(curLayer[i - 1][j][k]));
                      cell.addNeighbor(&(curLayer[i + 1][j][k]));

                      cell.addNeighbor(&(curLayer[i][j - 1][k]));
                      cell.addNeighbor(&(curLayer[i][j + 1][k]));

                      cell.addNeighbor(&(curLayer[i][j][k - 1]));
                      cell.addNeighbor(&(curLayer[i][j][k + 1]));
                  } else {
                      if (i == 0) {
                          cell.addNeighbor(&nullCell);
                          cell.addNeighbor(&(curLayer[i + 1][j][k]));
                      } else if (i == nx - 1) {
                          cell.addNeighbor(&(curLayer[i - 1][j][k]));
                          cell.addNeighbor(&nullCell);
                      } else {
                          cell.addNeighbor(&(curLayer[i - 1][j][k]));
                          cell.addNeighbor(&(curLayer[i + 1][j][k]));
                      }

                      if (j == 0) {
                          cell.addNeighbor(&nullCell);
                          cell.addNeighbor(&(curLayer[i][j + 1][k]));
                      } else if (j == ny - 1) {
                          cell.addNeighbor(&(curLayer[i][j - 1][k]));
                          cell.addNeighbor(&nullCell);
                      } else {
                          cell.addNeighbor(&(curLayer[i][j - 1][k]));
                          cell.addNeighbor(&(curLayer[i][j + 1][k]));
                      }

                      if (k == 0) {
                          cell.addNeighbor(&nullCell);
                          cell.addNeighbor(&(curLayer[i][j][k + 1]));
                      } else if (k == nz - 1) {
                          cell.addNeighbor(&(curLayer[i][j][k - 1]));
                          cell.addNeighbor(&nullCell);
                      } else {
                          Cell<T>& zero = curLayer[0][0][0];
                          cell.addNeighbor(&(curLayer[i][j][k - 1]));
                          cell.addNeighbor(&(curLayer[i][j][k + 1]));
                      }
                  }
              }
          }
      }
    return curLayer;
}

template <class T>
Cell<Point<T>>*** generateGasCells(
        long int nx, long int ny, long int nz,
        T hx, T hy, T hz,
        T tau, T t0, T t1,
        long int paramCount
) {
    Cell<Point<T>>*** curLayer = create3DArray<Cell<Point<T>>>(nx, ny, nz);

    int center = std::floor((double)nx/2);
    T E = 0;
    T e = 0;
    Point<T> leftFiller = Point<T>(paramCount);
    e = 1 / (0.4);
    leftFiller[0] = 1;
    leftFiller[1] = 0;
    leftFiller[2] = 0;
    leftFiller[3] = 0;
    leftFiller[4] = e + (pow(leftFiller[1], 2) + pow(leftFiller[2],2) + pow(leftFiller[3],2))/2;
    leftFiller[5] = 1;

    Point<T> rightFiller = Point<T>(paramCount);
    e = 0.1 / (0.4);
    rightFiller[0] = 0.125;
    rightFiller[1] = 0;
    rightFiller[2] = 0;
    rightFiller[3] = 0;
    rightFiller[4] = e + (pow(leftFiller[1], 2) + pow(leftFiller[2],2) + pow(leftFiller[3],2))/2;
    rightFiller[5] = 0.1;

    for (int i = 0; i < nx; i++) {
        if(i <= center) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    Cell<Point<T>> cell = Cell<Point<T>>(1, leftFiller, 3);
                    curLayer[i][j][k] = cell;
                }
            }
        } else {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    curLayer[i][j][k] = Cell<Point<T>>(1, rightFiller, 3);
                }
            }
        }
    }

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k].setType(0);
            curLayer[nx-1][j][k].setType(1);
        }
    }

    for(int i = 1; i < nx-1; i++) {
        for(int k = 0; k < nz; k++) {
            curLayer[i][0][k].setType(2);
            curLayer[i][ny-1][k].setType(2);
        }
    }
    for(int i = 1; i < nx - 1; i++) {
        for(int j = 0; j < ny; j++){
            curLayer[i][j][0].setType(2);
            curLayer[i][j][nz-1].setType(2);
        }
    }
    Cell<Point<T>>& nullCell = *(new Cell<Point<T>>(1, Point<T>(paramCount, 0.0), 1));
    return curLayer;
}

double*** convertToArray(Cell<double>*** cells, long int nx, long int ny, long int nz) {
    double*** result = create3DArray<double>(nx, ny, nz);
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                result[i][j][k] = cells[i][j][k][0];
            }
        }
    }
    return result;
}

double*** convertToArray(Cell<Point<double>>*** cells, long int nx, long int ny, long int nz) {
    double*** result = create3DArray<double>(nx, ny, nz);
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                result[i][j][k] = cells[i][j][k][0][1] / cells[i][j][k][0][0];
            }
        }
    }
    return result;
}


void printSodU(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][1] / cells[i][1][1][0][0]<<std::endl;
    }
    file.close();
}

void printSodV(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][2] / cells[i][1][1][0][0]<<std::endl;
    }
    file.close();
}

void printSodW(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][3] / cells[i][1][1][0][0]<<std::endl;
    }
    file.close();
}

void printSodE(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][4] / cells[i][1][1][0][0]<<std::endl;
    }
    file.close();
}

void printSodP(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][5]<<std::endl;
    }
    file.close();
}

void printSodRo(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][1][1][0][0]<<std::endl;
    }
    file.close();
}


int main() {
    long int Nx = 100;
    long int Ny = 5;
    long int Nz = 5;
    double x0 = 0; double x1 = 1;
    double y0 = 0; double y1 = 1;
    double z0 = 0; double z1 = 1;

    double hx = (x1 - x0) / Nx;
    double hy = (y1 - y0) / Ny;
    double hz = (z1 - z0) / Nz;

    double hmin = std::min(std::min(hx, hy), hz);

    double tau = pow(hmin,2);
    double t0 = 0; double t1 = 0.2;

    Cell<Point<double>>*** cells = generateGasCells(
            Nx, Ny, Nz,
            hx, hy, hz,
            tau, t0, t1, 6
    );
    cells = gasDynamic(cells,
               Nx, Ny, Nz,
               x0, x1,
               y0, y1,
               z0, z1,
               hx, hy, hz,
               tau, t0, t1,
               1,
               100,
               "Sod");
    double*** result = convertToArray(cells, Nx, Ny, Nz);


    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
    printSodU(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "U.txt");
    printSodV(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "V.txt");
    printSodW(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "W.txt");
    printSodE(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "E.txt");
    printSodP(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "P.txt");
    printSodRo(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "Ro.txt");

    std::cout << "RESULT SOLVE" << std::endl;
    return 0;
}