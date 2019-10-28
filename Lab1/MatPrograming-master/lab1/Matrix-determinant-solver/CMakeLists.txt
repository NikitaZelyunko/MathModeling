cmake_minimum_required(VERSION 3.0.0)

project(MatrixDeterminantSolver)

set(INCLUDE_DIR include)
set(SOURCE_DIR src) 

add_library(MatrixDeterminantSolver STATIC 
    ${INCLUDE_DIR}/MatrixDeterminantSolver.h
    ${SOURCE_DIR}/MatrixDeterminantSolver.cpp)

target_link_libraries(MatrixDeterminantSolver Utils)
target_link_libraries(MatrixDeterminantSolver MatrixAbstractSolver)
target_link_libraries(MatrixDeterminantSolver Matrix)

target_include_directories(MatrixDeterminantSolver PRIVATE ${SOURCE_DIR})
target_include_directories(MatrixDeterminantSolver PUBLIC ${INCLUDE_DIR})