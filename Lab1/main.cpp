
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


int main() {
    long int N = 30;
    double x0 = 0; double x1 = 1;
    double y0 = 0; double y1 = 1;
    double z0 = 0; double z1 = 1;

    double tau = 0.001; double t0 = 0; double t1 = 2;

//    double ***result = Reshenie_Uravn_Teploprovodnosti_methodom_progonki_yavn<double>(
//            N,
//            x0, x1,
//            y0, y1,
//            z0, z1,
//            tau, t0, t1
//    );

    double ***result = Reshenie_Uravn_Teploprovodnosti_flux<double>(
            N,
            x0, x1,
            y0, y1,
            z0, z1,
            tau, t0, t1
    );

    std::cout<<"RESULT SOLVE"<<std::endl;
    printCSV(result, N, x0, x1, y0, y1, z0, z1, "resultFlux.csv");
    delete3DArray(result, N, N);
    return 0;
}