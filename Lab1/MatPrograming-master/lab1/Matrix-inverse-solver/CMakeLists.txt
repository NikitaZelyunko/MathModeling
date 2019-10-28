cmake_minimum_required(VERSION 3.0.0)

project(MatrixInverseSolver)

set(INCLUDE_DIR include)
set(SOURCE_DIR src)

add_library(MatrixInverseSolver STATIC 
    ${INCLUDE_DIR}/MatrixInverseSolver.h
    ${SOURCE_DIR}/MatrixInverseSolver.cpp
)

target_link_libraries(MatrixInverseSolver Utils)
target_link_libraries(MatrixInverseSolver Matrix)
target_link_libraries(MatrixInverseSolver MatrixAbstractSolver)
target_link_libraries(MatrixInverseSolver MatrixDeterminantSolver)

target_include_directories(MatrixInverseSolver PRIVATE ${SOURCE_DIR})
target_include_directories(MatrixInverseSolver PUBLIC ${INCLUDE_DIR})