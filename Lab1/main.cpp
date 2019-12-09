
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
Cell<Point<T>>*** generateSodCells(
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
//    solveNeighbors<Point<T>>(curLayer, nullCell, nx, ny, nz, paramCount);
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

template<class T>
bool isInSphere(T x, T y, T z, T x0, T y0, T z0, T r) {
    return (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) <= pow(r, 2);
}

template <class T>
Cell<Point<T>>*** generateHeliumCylinderCells(
        T hx, T hy, T hz,
        long int nx, long int ny, long int nz,
        T tau, T t0, T t1,
        T D,
        long int paramCount
) {
    T R = 0.5 * D;

    T x0 = 0, x1 = 6.5 * D;
    T y0 = 0, y1 = 1.78 * D;
    T z0 = 0, z1 = 1.78 * D;

    T postShockX0 = x0, postShockX1 = 2 * D;
    T postShockY0 = y0, postShockY1 = y1;
    T postShockZ0 = z0, postShockZ1 = z1;

    T sphereCenterX = 3.5 * D, sphereCenterY = y1 / 2, sphereCenterZ = z1 / 2;

    T preShockFirstX0 = postShockX1, preShockFirstX1 = preShockFirstX0 + 0.5 * D;
    T preShockFirstY0 = y0, preShockFirstY1 = y1;
    T preShockFirstZ0 = z0, preShockFirstZ1 = z1;

    T preShockSecondX0 = preShockFirstX1, preShockSecondX1 = preShockSecondX0 + D;
    T preShockSecondY0 = y0, preShockSecondY1 = y1 / 2 - R;
    T preShockSecondZ0 = z0, preShockSecondZ1 = z1 / 2 - R;

    T preShockThirdX0 = preShockFirstX1, preShockThirdX1 = preShockThirdX0 + D;
    T preShockThirdY0 = y1 / 2 + R, preShockThirdY1 = y1;
    T preShockThirdZ0 = z1 / 2 + R, preShockThirdZ1 = z1;

    T preShockFourthX0 = preShockThirdX1, preShockFourthX1 = x1;
    T preShockFourthY0 = y0, preShockFourthY1 = y1;
    T preShockFourthZ0 = z0, preShockFourthZ1 = z1;

    T bubbleSquareX0 = preShockFirstX1, bubbleSquareX1 = preShockFirstX1 + D;
    T bubbleSquareY0 = preShockSecondY1, bubbleSquareY1 = bubbleSquareY0 + R;
    T bubbleSquareZ0 = preShockSecondZ1, bubbleSquareZ1 = bubbleSquareZ0 + R;


    int firstSectorXLength = floor((postShockX1 - postShockX0) / hx), firstSectorYLength = floor((postShockY1 - postShockY0)/ hy), firstSectorZLength = floor((postShockZ1 - postShockZ0) / hz);
    int
    firstSectorX0 = 0, firstSectorX1 = firstSectorX0 + firstSectorXLength,
    firstSectorY0 = 0, firstSectorY1 = ny,
    firstSectorZ0 = 0, firstSectorZ1 = nz;

    int secondSectorXLength = floor((preShockFirstX1 - preShockFirstX0) /hx), secondSectorYLength = floor((preShockFirstY1 - preShockFirstY0) /hy), secondSectorZLength = floor((preShockFirstZ1 - preShockFirstZ0) /hz);
    int
    secondSectorX0 = firstSectorX1, secondSectorX1 = secondSectorX0 + secondSectorXLength,
    secondSectorY0 = 0, secondSectorY1 = ny,
    secondSectorZ0 = 0, secondSectorZ1 = nz;

    int thirdSectorXLength = floor((preShockSecondX1 - preShockSecondX0) / hx), thirdSectorYLength = floor((preShockSecondY1 - preShockSecondY0) / hy), thirdSectorZLength = floor((preShockSecondZ1 - preShockSecondZ0) / hz);
    int
    thirdSectorX0 = secondSectorX1, thirdSectorX1 = thirdSectorX0 + thirdSectorXLength,
    thirdSectorY0 = 0, thirdSectorY1 = thirdSectorY0 + thirdSectorYLength,
    thirdSectorZ0 = 0, thirdSectorZ1 = thirdSectorZ0 + thirdSectorZLength;

    int bubbleSectorXLength = floor((bubbleSquareX1 - bubbleSquareX0) / hx), bubbleSectorYLength = floor((bubbleSquareY1 - bubbleSquareY0) / hy), bubbleSectorZLength = floor((bubbleSquareZ1 - bubbleSquareZ0) / hz);
    int
    bubbleSectorX0 = secondSectorX1, bubbleSectorX1 = bubbleSectorX0 + bubbleSectorXLength,
    bubbleSectorY0 = thirdSectorY1, bubbleSectorY1 = bubbleSectorY0 + bubbleSectorYLength,
    bubbleSectorZ0 = thirdSectorZ1, bubbleSectorZ1 = bubbleSectorZ0 + bubbleSectorZLength;

    int fourthSectorXLength = floor((preShockThirdX1 - preShockThirdX0) / hx), fourthSectorYLength = floor((preShockThirdY1 - preShockThirdY0) / hy), fourthSectorZLength = floor((preShockThirdZ1 - preShockThirdZ0) / hz);
    int
    fourthSectorX0 = secondSectorX1, fourthSectorX1 = fourthSectorX0 + fourthSectorXLength,
    fourthSectorY0 = bubbleSectorY1, fourthSectorY1 = ny,
    fourthSectorZ0 = bubbleSectorZ1, fourthSectorZ1 = nz;

    int fifthSectorXLength = floor((preShockFourthX1 - preShockFourthX0) / hx), fifthSectorYLength = floor((preShockFourthY1 - preShockFourthY0) / hy), fifthSectorZLength = floor((preShockFourthZ1 - preShockFourthZ0) / hz);
    int
    fifthSectorX0 = fourthSectorX1, fifthSectorX1 = nx,
    fifthSectorY0 = fourthSectorY1, fifthSectorY1 = ny,
    fifthSectorZ0 = fourthSectorZ1, fifthSectorZ1 = nz;

    Cell<Point<T>>*** curLayer = create3DArray<Cell<Point<T>>>(nx, ny, nz);

    T E = 0;
    T e = 0;
    Cell<T> postShockAir = Cell<T>(paramCount);
    e = 1 / (0.4);
    postShockAir[0] = 1.65;
    postShockAir[1] = 114.4 * postShockAir[0];
    postShockAir[2] = 0;
    postShockAir[3] = 0;
    postShockAir[4] = e + (pow(postShockAir[1], 2) + pow(postShockAir[2],2) + pow(postShockAir[3],2))/2;
    postShockAir[5] = 158900;
    postShockAir.setType(3);


    Cell<T> preShockAir = Cell<T>(paramCount);
    e = 0.1 / (0.4);
    preShockAir[0] = 1.20;
    preShockAir[1] = 0;
    preShockAir[2] = 0;
    preShockAir[3] = 0;
    preShockAir[4] = e + (pow(preShockAir[1], 2) + pow(preShockAir[2],2) + pow(preShockAir[3],2))/2;
    preShockAir[5] = 101325;
    preShockAir.setType(3);

    Cell<T> helium = Cell<T>(paramCount);
    e = 0.1 / (0.4);
    helium[0] = 0.166;
    helium[1] = 0;
    helium[2] = 0;
    helium[3] = 0;
    helium[4] = e + (pow(helium[1], 2) + pow(helium[2],2) + pow(helium[3],2))/2;
    helium[5] = 101325;
    helium.setType(3);

    for(int i = firstSectorX0; i < firstSectorX1; i++) {
        for(int j = firstSectorY0; j < firstSectorY1; j++) {
            for(int k = firstSectorZ0; k < firstSectorZ1; k++) {
                curLayer[i][j][k] = postShockAir;
            }
        }
    }

    for(int i = secondSectorX0; i < secondSectorX1; i++) {
        for(int j = secondSectorY0; j < secondSectorY1; j++) {
            for(int k = secondSectorZ0; k < secondSectorZ1; k++) {
                curLayer[i][j][k] = preShockAir;
            }
        }
    }

    for(int i = thirdSectorX0; i < thirdSectorX1; i++) {
        for(int j = thirdSectorY0; j < thirdSectorY1; j++) {
            for(int k = thirdSectorZ0; k < thirdSectorZ1; k++) {
                curLayer[i][j][k] = preShockAir;
            }
        }
    }

    for(int i = bubbleSectorX0; i < bubbleSectorX1; i++) {
        for(int j = bubbleSectorY0; j < bubbleSectorY1; j++) {
            for(int k = bubbleSectorZ0; k < bubbleSectorZ1; k++) {

                T
                xCellLeft = x0 + i * hx, xCellRight = xCellLeft + hx,
                yCellLeft = y0 + j * hy, yCellRight = yCellLeft + hy,
                zCellLeft = z0 + k * hz, zCellRight = zCellLeft + hz;

                if(isInSphere(xCellLeft + (xCellRight - xCellLeft) / 2, yCellLeft + (yCellRight - yCellLeft) / 2, zCellLeft + (zCellRight - zCellLeft) / 2, sphereCenterX, sphereCenterY, sphereCenterZ, R)) {
                    curLayer[i][j][k] = helium;
                } else {
                    curLayer[i][j][k] = preShockAir;
                }
            }
        }
    }

    for(int i = fourthSectorX0; i < fourthSectorX1; i++) {
        for(int j = fourthSectorY0; j < fourthSectorY1; j++) {
            for(int k = fourthSectorZ0; k < fourthSectorZ1; k++) {
                curLayer[i][j][k] = preShockAir;
            }
        }
    }

    for(int i = fifthSectorX0; i < fifthSectorX1; i++) {
        for(int j = fifthSectorY0; j < fifthSectorY1; j++) {
            for(int k = fifthSectorZ0; k < fifthSectorZ1; k++) {
                curLayer[i][j][k] = preShockAir;
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
    return curLayer;
}

double*** convertRoToArray(Cell<Point<double>>*** cells, long int nx, long int ny, long int nz) {
    double*** result = create3DArray<double>(nx, ny, nz);
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                result[i][j][k] = cells[i][j][k][0][0];
            }
        }
    }
    return result;
}


int main() {

    double D = 0.05;
    long int Nx = 100;
    long int Ny = 30;
    long int Nz = 30;
    double x0 = 0; double x1 = 6.5 * D;
    double y0 = 0; double y1 = 1.78 * D;
    double z0 = 0; double z1 = 1.78 * D;


    double hx = (x1 - x0) / Nx;
    double hy = (y1 - y0) / Ny;
    double hz = (z1 - z0) / Nz;

    double hmin = std::min(std::min(hx, hy), hz);

    double tau = pow(hmin,2);
//    double tau = 2.5 / 100000;
    double t0 = 0; double t1 = 0.02;

//    double ***result = Reshenie_Uravn_Teploprovodnosti_methodom_progonki_yavn<double>(
//        Nx,
//        x0, x1,
//        y0, y1,
//        z0, z1,
//        tau, t0, t1
//    );
//
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, "Result1.vts");
//
//    double ***result2 = Reshenie_Uravn_Teploprovodnosti_flux<double>(
//            Nx,
//            x0, x1,
//            y0, y1,
//            z0, z1,
//            tau, t0, t1
//    );
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, "Result2.vts");

//    Cell<double>*** cells = generateThermalCells(
//                Nx, Ny, Nz,
//                hx, hy, hz,
//                tau, t0, t1, 1
//            );
//    std::cout<<"All right"<<std::endl;
//    gasDynamic(cells,
//            Nx, Ny, Nz,
//            hx, hy, hz,
//            tau, t0, t1,
//            1);
//    std::cout<<"All right"<<std::endl;
//    double*** result = convertToArray(cells, Nx, Ny, Nz);
//
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, "TestGasDynamic.vts");
//
//    std::cout << "RESULT SOLVE" << std::endl;

    Cell<Point<double>>*** cells = generateHeliumCylinderCells(
            hx, hy, hz,
            Nx, Ny, Nz,
            tau, t0, t1,
            D,
            6
    );
//    std::cout<<"All right"<<std::endl;
//    cells = gasDynamic(cells,
//               Nx, Ny, Nz,
//               hx, hy, hz,
//               tau, t0, t1,
//               1);
//    std::cout<<"All right"<<std::endl;
//    for(int i=0; i < Nx; i++) {
//        for(int j = 0; j < Ny; j++) {
//            for(int k = 0; k < Nz; k++) {
//                std::cout<<i<<","<<j<<","<<k;
//                cells[i][j][k][0].print("u:");
//            }
//        }
//    }
    double*** result = convertRoToArray(cells, Nx, Ny, Nz);
    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestHeliumBubble.vts");

//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    printSodU(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "U.txt");
//    printSodV(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "V.txt");
//    printSodW(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "W.txt");
//    printSodE(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "E.txt");
//    printSodP(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "P.txt");
//    printSodRo(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "Ro.txt");

    std::cout << "RESULT SOLVE" << std::endl;
//    double ***tochn = create3DArray<double>(N, N, N);
//    double ***razn = create3DArray<double>(N, N, N);
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            for (int k = 0; k < N; k++) {
//                tochn[i][j][k] = pow(2.7, -3 * pow(3.14, 2) * t1) * sin(3.14 * (x0 + i*h)) * sin(3.14 * (y0 + j*h)) * sin(3.14 * (z0 + k*h));
//
//
//                razn[i][j][k] = 0;
//            }
//        }
//    }
//    for (int i = 1; i < N - 1; i++) {
//        for (int j = 1; j < N - 1; j++) {
//            for (int k = 1; k < N - 1; k++) {
//                razn[i][j][k] = abs(result[i][j][k] - tochn[i][j][k]);
//                if (razn[i][j][k] > 0.1) {
//                    std::cout << i << "," << j << "," << k << " " << result[i][j][k] << " " << tochn[i][j][k] << std::endl;
//                }
//            }
//        }
//    }
//    printCSV(result, N, x0, x1, y0, y1, z0, z1, "resultFlux.csv");
//    printCSV(tochn, N, x0, x1, y0, y1, z0, z1, "tochnFlux.csv");
//    printCSV(razn, N, x0, x1, y0, y1, z0, z1, "raznFlux.csv");
//    delete3DArray(result, Nx, Ny);
//    delete3DArray(tochn, N, N);
//    delete3DArray(razn, N, N);
    return 0;
}