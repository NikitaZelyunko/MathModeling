cmake_minimum_required(VERSION 3.14)
project(Lab1)

set(CMAKE_CXX_STANDARD 14)

add_executable(Lab1 main.cpp header.h VTSFormater.h Cell.cpp Cell.h)

add_subdirectory(MatPrograming-master/lab1)

target_link_libraries(Lab1 MatProgramming)