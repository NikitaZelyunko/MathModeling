
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


int main() {
    long int Nx = 20;
    long int Ny = 20;
    long int Nz = 20;
    double x0 = 0; double x1 = 1;
    double y0 = 0; double y1 = 1;
    double z0 = 0; double z1 = 1;

    double hx = (x1 - x0) / Nx;
    double hy = (y1 - y0) / Ny;
    double hz = (z1 - z0) / Nz;

    double hmin = std::min(std::min(hx, hy), hz);

    double tau = pow(hmin,2)/10/3; double t0 = 0; double t1 = 0.1;

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

    Cell<double>*** cells = generateThermalCells(
                Nx, Ny, Nz,
                hx, hy, hz,
                tau, t0, t1, 1
            );
    std::cout<<"All right"<<std::endl;
    gasDynamic(cells,
            Nx, Ny, Nz,
            hx, hy, hz,
            tau, t0, t1,
            1);
    std::cout<<"All right"<<std::endl;
    double*** result = convertToArray(cells, Nx, Ny, Nz);

    VTSFormateer(result, Nx, Ny, Nz, x0, x1, y0, y1, z0, z1, "TestGasDynamic.vts");

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
    delete3DArray(result, Nx, Ny);
//    delete3DArray(tochn, N, N);
//    delete3DArray(razn, N, N);
    return 0;
}