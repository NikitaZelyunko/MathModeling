cmake_minimum_required(VERSION 3.0.0)

project(MatrixAbstractSolver)

set(INCLUDE_DIR include)
set(SOURCE_DIR src)

add_library(MatrixAbstractSolver STATIC
    ${INCLUDE_DIR}/MatrixAbstractSolver.h
    ${SOURCE_DIR}/MatrixAbstractSolver.cpp
)

target_link_libraries(MatrixAbstractSolver Utils)
target_link_libraries(MatrixAbstractSolver AbstractSolver)
target_link_libraries(MatrixAbstractSolver Matrix)

target_include_directories(MatrixAbstractSolver PRIVATE ${SOURCE_DIR})
target_include_directories(MatrixAbstractSolver PUBLIC ${INCLUDE_DIR})