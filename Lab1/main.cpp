
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
    rightFiller[4] = e + (pow(rightFiller[1], 2) + pow(rightFiller[2],2) + pow(rightFiller[3],2))/2;
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

    Cell<Point<T>>* leftNeighbour = new Cell<Point<T>>(1, leftFiller, 0);
    Cell<Point<T>>* rightNeighbour = new Cell<Point<T>>(1, rightFiller, 1);

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k].setType(0);
            curLayer[0][j][k].addNeighbor(leftNeighbour);
            curLayer[nx-1][j][k].setType(1);
            curLayer[nx-1][j][k].addNeighbor(rightNeighbour);
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
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][1] / cells[i][0][0][0][0]<<std::endl;
    }
    file.close();
}

void printSodV(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][2] / cells[i][0][0][0][0]<<std::endl;
    }
    file.close();
}

void printSodW(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][3] / cells[i][0][0][0][0]<<std::endl;
    }
    file.close();
}

void printSodE(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][4] / cells[i][0][0][0][0]<<std::endl;
    }
    file.close();
}

void printSodP(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][5]<<std::endl;
    }
    file.close();
}

void printSodRo(Cell<Point<double>>*** cells, int nx, int ny, int nz, double x0, double y0, double z0, double hx, double hy, double hz, std::string filename) {
    std::ofstream file(filename);
    for(int i = 0; i < nx; i++) {
        file<<std::fixed<<std::setw(6)<<std::setprecision(5)<<x0 + hx*i<<" "<<cells[i][0][0][0][0]<<std::endl;
    }
    file.close();
}

template<class T>
bool isInSphere(T x, T y, T z, T x0, T y0, T z0, T r) {
    return pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2) <= pow(r, 2);
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

    T sphereCenterX = 3 * D, sphereCenterY = y1 / 2, sphereCenterZ = z1 / 2;

    T preShockFirstX0 = postShockX1, preShockFirstX1 = preShockFirstX0 + 0.5 * D;
    T preShockFirstY0 = y0, preShockFirstY1 = y1;
    T preShockFirstZ0 = z0, preShockFirstZ1 = z1;

    T bubbleNeighbourBottomX0 = preShockFirstX1, bubbleNeighbourBottomX1 = bubbleNeighbourBottomX0 + D;
    T bubbleNeighbourBottomY0 = y0, bubbleNeighbourBottomY1 = y1 / 2 - R;
    T bubbleNeighbourBottomZ0 = z0, bubbleNeighbourBottomZ1 = z1;

    T bubbleSquareX0 = preShockFirstX1, bubbleSquareX1 = bubbleSquareX0 + D;
    T bubbleSquareY0 = bubbleNeighbourBottomY1, bubbleSquareY1 = bubbleSquareY0 + D;
    T bubbleSquareZ0 = z1 / 2 - R, bubbleSquareZ1 = bubbleSquareZ0 + D;

    T bubbleNeighbourTopX0 = preShockFirstX1, bubbleNeighbourTopX1 = bubbleNeighbourTopX0 + D;
    T bubbleNeighbourTopY0 = bubbleSquareY1, bubbleNeighbourTopY1 = y1;
    T bubbleNeighbourTopZ0 = z0, bubbleNeighbourTopZ1 = z1;

    T bubbleNeighbourFrontX0 = preShockFirstX1, bubbleNeighbourFrontX1 = bubbleNeighbourFrontX0 + D;
    T bubbleNeighbourFrontY0 = bubbleSquareY0, bubbleNeighbourFrontY1 = bubbleSquareY1;
    T bubbleNeighbourFrontZ0 = z0, bubbleNeighbourFrontZ1 = bubbleSquareZ0;

    T bubbleNeighbourBackX0 = preShockFirstX1, bubbleNeighbourBackX1 = bubbleNeighbourBackX0 + D;
    T bubbleNeighbourBackY0 = bubbleSquareY0, bubbleNeighbourBackY1 = bubbleSquareY1;
    T bubbleNeighbourBackZ0 = bubbleSquareZ1, bubbleNeighbourBackZ1 = z1;

    T preShockThirdX0 = bubbleSquareX1, preShockThirdX1 = x1;
    T preShockThirdY0 = y0, preShockThirdY1 = y1;
    T preShockThirdZ0 = z0, preShockThirdZ1 = z1;


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

    int
            bubbleNeighbourBottomXLength = floor((bubbleNeighbourBottomX1 - bubbleNeighbourBottomX0) / hx) + 1,
            bubbleNeighbourBottomYLength = floor((bubbleNeighbourBottomY1 - bubbleNeighbourBottomY0) / hy) + 1,
            bubbleNeighbourBottomZLength = floor((bubbleNeighbourBottomZ1 - bubbleNeighbourBottomZ0) / hz) + 1;

    int
            bubbleNeighbourBottomSectorX0 = secondSectorX1, bubbleNeighbourBottomSectorX1 = bubbleNeighbourBottomSectorX0 + bubbleNeighbourBottomXLength,
            bubbleNeighbourBottomSectorY0 = 0, bubbleNeighbourBottomSectorY1 = bubbleNeighbourBottomSectorY0 + bubbleNeighbourBottomYLength,
            bubbleNeighbourBottomSectorZ0 = 0, bubbleNeighbourBottomSectorZ1 = nz;

    int
            bubbleNeighbourFrontXLength = floor((bubbleNeighbourFrontX1 - bubbleNeighbourFrontX0) / hx) + 1,
            bubbleNeighbourFrontYLength = floor((bubbleNeighbourFrontY1 - bubbleNeighbourFrontY0) / hy) + 1,
            bubbleNeighbourFrontZLength = floor((bubbleNeighbourFrontZ1 - bubbleNeighbourFrontZ0) / hz) + 1;

    int
            bubbleNeighbourFrontSectorX0 = secondSectorX1, bubbleNeighbourFrontSectorX1 = bubbleNeighbourFrontSectorX0 + bubbleNeighbourFrontXLength,
            bubbleNeighbourFrontSectorY0 = bubbleNeighbourBottomSectorY1, bubbleNeighbourFrontSectorY1 = bubbleNeighbourFrontSectorY0 + bubbleNeighbourFrontYLength,
            bubbleNeighbourFrontSectorZ0 = 0, bubbleNeighbourFrontSectorZ1 = bubbleNeighbourFrontSectorZ0 + bubbleNeighbourFrontZLength;

    int
            bubbleNeighbourTopXLength = floor((bubbleNeighbourTopX1 - bubbleNeighbourTopX0) / hx) + 1,
            bubbleNeighbourTopYLength = floor((bubbleNeighbourTopY1 - bubbleNeighbourTopY0) / hy) + 1,
            bubbleNeighbourTopZLength = floor((bubbleNeighbourTopZ1 - bubbleNeighbourTopZ0) / hz) + 1;

    int
            bubbleNeighbourTopSectorX0 = secondSectorX1, bubbleNeighbourTopSectorX1 = bubbleNeighbourTopSectorX0 + bubbleNeighbourTopXLength,
            bubbleNeighbourTopSectorY0 = bubbleNeighbourFrontSectorY1, bubbleNeighbourTopSectorY1 = ny,
            bubbleNeighbourTopSectorZ0 = 0, bubbleNeighbourTopSectorZ1 = nz;

    int
            bubbleSectorXLength = floor((bubbleSquareX1 - bubbleSquareX0) / hx) + 1,
            bubbleSectorYLength = floor((bubbleSquareY1 - bubbleSquareY0) / hy) + 1,
            bubbleSectorZLength = floor((bubbleSquareZ1 - bubbleSquareZ0) / hz) + 1;
    int
            bubbleSectorX0 = secondSectorX1, bubbleSectorX1 = bubbleSectorX0 + bubbleSectorXLength,
            bubbleSectorY0 = bubbleNeighbourBottomSectorY1, bubbleSectorY1 = bubbleSectorY0 + bubbleSectorYLength,
            bubbleSectorZ0 = bubbleNeighbourFrontSectorZ1, bubbleSectorZ1 = bubbleSectorZ0 + bubbleSectorZLength;

    int
            bubbleNeighbourBackXLength = floor((bubbleNeighbourBackX1 - bubbleNeighbourBackX0) / hx) + 1,
            bubbleNeighbourBackYLength = floor((bubbleNeighbourBackY1 - bubbleNeighbourBackY0) / hy) + 1,
            bubbleNeighbourBackZLength = floor((bubbleNeighbourBackZ1 - bubbleNeighbourBackZ0) / hz) + 1;

    int
            bubbleNeighbourBackSectorX0 = secondSectorX1, bubbleNeighbourBackSectorX1 = bubbleNeighbourBackSectorX0 + bubbleNeighbourBackXLength,
            bubbleNeighbourBackSectorY0 = bubbleNeighbourBottomSectorY1, bubbleNeighbourBackSectorY1 = bubbleNeighbourBackSectorY0 + bubbleNeighbourBackYLength,
            bubbleNeighbourBackSectorZ0 = bubbleSectorZ1, bubbleNeighbourBackSectorZ1 = nz;

    int
            thirdSectorXLength = floor((preShockThirdX1 - preShockThirdX0) / hx),
            thirdSectorYLength = floor((preShockThirdY1 - preShockThirdY0) / hy),
            thirdSectorZLength = floor((preShockThirdZ1 - preShockThirdZ0) / hz);

    int
            thirdSectorX0 = bubbleNeighbourBackSectorX1, thirdSectorX1 = nx,
            thirdSectorY0 = 0, thirdSectorY1 = ny,
            thirdSectorZ0 = 0, thirdSectorZ1 = nz;

    Cell<Point<T>>*** curLayer = create3DArray<Cell<Point<T>>>(nx, ny, nz);

    T
            ro = 0,
            p = 0,
            E = 0,
            e = 0,
            sigma = 0,
            u = 0, v = 0, w = 0;
    Point<T> postShockAir = Point<T>(paramCount);


    ro = 1.65;
    p = 158900;
    sigma = 1.4;
    u = 114.4, v = 0, w = 0;
    e = p / (ro * (sigma - 1));
    E = e + (pow(u, 2) + pow(v, 2) + pow(w, 2)) / 2;

    postShockAir[0] = ro;
    postShockAir[1] = ro * u;
    postShockAir[2] = ro * v;
    postShockAir[3] = ro * w;
    postShockAir[4] = ro * E;
    postShockAir[5] = p;

    ro = 1.20;
    p = 101325;
    sigma = 1.4;
    u = 0, v = 0, w = 0;
    e = p / (ro * (sigma - 1));
    E = e + (pow(u, 2) + pow(v, 2) + pow(w, 2)) / 2;

    Point<T> preShockAir = Point<T>(paramCount);

    preShockAir[0] = ro;
    preShockAir[1] = ro * u;
    preShockAir[2] = ro * v;
    preShockAir[3] = ro * w;
    preShockAir[4] = ro * E;
    preShockAir[5] = p;

    ro = 0.166;
    p = 101325;
    sigma = 1.4;
    u = 0, v = 0, w = 0;
    e = p / (ro * (sigma - 1));
    E = e + (pow(u, 2) + pow(v, 2) + pow(w, 2)) / 2;

    Point<T> helium = Point<T>(paramCount);

    helium[0] = ro;
    helium[1] = ro * u;
    helium[2] = ro * v;
    helium[3] = ro * w;
    helium[4] = ro * E;
    helium[5] = p;


    for(int i = firstSectorX0; i < firstSectorX1; i++) {
        for(int j = firstSectorY0; j < firstSectorY1; j++) {
            for(int k = firstSectorZ0; k < firstSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, postShockAir, 3);
            }
        }
    }

    for(int i = secondSectorX0; i < secondSectorX1; i++) {
        for(int j = secondSectorY0; j < secondSectorY1; j++) {
            for(int k = secondSectorZ0; k < secondSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
            }
        }
    }

    for(int i = bubbleNeighbourBottomSectorX0; i < bubbleNeighbourBottomSectorX1; i++) {
        for(int j = bubbleNeighbourBottomSectorY0; j < bubbleNeighbourBottomSectorY1; j++) {
            for(int k = bubbleNeighbourBottomSectorZ0; k < bubbleNeighbourBottomSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
            }
        }
    }

    for(int i = bubbleNeighbourFrontSectorX0; i < bubbleNeighbourFrontSectorX1; i++) {
        for(int j = bubbleNeighbourFrontSectorY0; j < bubbleNeighbourFrontSectorY1; j++) {
            for(int k = bubbleNeighbourFrontSectorZ0; k < bubbleNeighbourFrontSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
            }
        }
    }

    for(int i = bubbleNeighbourTopSectorX0; i < bubbleNeighbourTopSectorX1; i++) {
        for(int j = bubbleNeighbourTopSectorY0; j < bubbleNeighbourTopSectorY1; j++) {
            for(int k = bubbleNeighbourTopSectorZ0; k < bubbleNeighbourTopSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
            }
        }
    }

    for(int i = bubbleNeighbourBackSectorX0; i < bubbleNeighbourBackSectorX1; i++) {
        for(int j = bubbleNeighbourBackSectorY0; j < bubbleNeighbourBackSectorY1; j++) {
            for(int k = bubbleNeighbourBackSectorZ0; k < bubbleNeighbourBackSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
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
                    curLayer[i][j][k] = Cell<Point<T>>(1, helium, 3);
                } else {
                    curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
                }
            }
        }
    }

    for(int i = thirdSectorX0; i < thirdSectorX1; i++) {
        for(int j = thirdSectorY0; j < thirdSectorY1; j++) {
            for(int k = thirdSectorZ0; k < thirdSectorZ1; k++) {
                curLayer[i][j][k] = Cell<Point<T>>(1, preShockAir, 3);
            }
        }
    }

    Cell<Point<T>>* postShockNeighbour = new Cell<Point<T>>(1, postShockAir, 0);

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            curLayer[0][j][k].setType(0);
            curLayer[0][j][k].addNeighbor(postShockNeighbour);
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


int main() {
//    long int Nx = 200;
//    long int Ny = 3;
//    long int Nz = 3;
//    double x0 = 0; double x1 = 1;
//    double y0 = 0; double y1 = 1;
//    double z0 = 0; double z1 = 1;
//
//    double hx = (x1 - x0) / Nx;
//    double hy = (y1 - y0) / Ny;
//    double hz = (z1 - z0) / Nz;
//
//    double hmin = std::min(std::min(hx, hy), hz);
//
//    double tau = pow(hmin,2);
//    double t0 = 0; double t1 = 0.2;
//
//    Cell<Point<double>>*** cells = generateGasCells(
//            Nx, Ny, Nz,
//            hx, hy, hz,
//            tau, t0, t1, 6
//    );
    double D = 0.05;
    long int Nx = 200;
    long int Ny = Nx / 6.5 * 1.78;
    long int Nz = Nx / 6.5 * 1.78;
    double x0 = 0; double x1 = 6.5 * D;
    double y0 = 0; double y1 = 1.78 * D;
    double z0 = 0; double z1 = 1.78 * D;

    double hx = (x1 - x0) / Nx;
    double hy = (y1 - y0) / Ny;
    double hz = (z1 - z0) / Nz;

    double hmin = std::min(std::min(hx, hy), hz);

    double tau = pow(hmin,2)/10;
    double t0 = 0; double t1 = 0.02;

    Cell<Point<double>>*** cells = generateHeliumCylinderCells(
            hx, hy, hz,
            Nx, Ny, Nz,
            tau, t0, t1,
            D,
            6
    );
    Cell<Point<double>>*** result = gasDynamic(cells,
                                               Nx, Ny, Nz,
                                               x0, x1,
                                               y0, y1,
                                               z0, z1,
                                               hx, hy, hz,
                                               tau, t0, t1,
                                               1,
                                               10,
                                               "HeliumBubble"
    );
    delete3DArray(cells, Nx, Ny);
//    cells = gasDynamic(cells,
//               Nx, Ny, Nz,
//               hx, hy, hz,
//               tau, t0, t1,
//               1);
//    double*** result = convertToArray(cells, Nx, Ny, Nz);
//
//
//    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, 6, 4, "TestGasDynamic.vts");
//    printSodU(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "U.txt");
//    printSodV(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "V.txt");
//    printSodW(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "W.txt");
//    printSodE(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "E.txt");
//    printSodP(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "P.txt");
//    printSodRo(cells, Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, "Ro.txt");

    std::cout << "RESULT SOLVE" << std::endl;
    return 0;
}