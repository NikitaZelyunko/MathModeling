#ifndef LAB1_VTSFORMATER_H
#define LAB1_VTSFORMATER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

template <class T>
void VTSFormateer(
        T*** arr,
        long int nx, long int ny, long int nz,
        T x0, T x1,
        T y0, T y1,
        T z0, T z1,
        int numberCount, int precision,
        std::string filename) {

    T hx = (x1 - x0) / nx;
    T hy = (y1 - y0) / ny;
    T hz = (z1 - z0) / nz;

    std::ofstream file(filename);

    file<<"<?xml version=\"1.0\"?>"<<std::endl;
    file<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;

    std::stringstream pieces;
    pieces<<"\"0 "<<nx-1<<" 0 "<<ny - 1<<" 0 "<<nz-1<<"\"";

    file<<"<StructuredGrid WholeExtent="<<pieces.str()<<">"<<std::endl;

    file<<"<Piece Extent="<<pieces.str()<<">"<<std::endl;
    file<<"<Points>"<<std::endl;
    file<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<std::endl;

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                file<<x0 + i*hx<<" "<<y0 + j*hy<<" "<<z0 + k*hz<<" ";
            }
            file<<std::endl;
        }
    }
    file<<"</DataArray>"<<std::endl;
    file<<"</Points>"<<std::endl;
    file<<"<CellData>"<<std::endl;
    file<<"</CellData>"<<std::endl;
    file<<"<PointData>"<<std::endl;

    file<<"<DataArray type=\"Float64\" Name=\"Temperature\" format=\"ascii\">"<<std::endl;
    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                file<<std::fixed<<std::setw(numberCount)<<std::setprecision(precision)<<arr[i][j][k]<<" ";
            }
            file<<std::endl;
        }
    }
    file<<"</DataArray>"<<std::endl;
    file<<"</PointData>"<<std::endl;
    file<<"</Piece>"<<std::endl;
    file<<"</StructuredGrid>"<<std::endl;
    file<<"</VTKFile>"<<std::endl;
}

#endif //LAB1_VTSFORMATER_H