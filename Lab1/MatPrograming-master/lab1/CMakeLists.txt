cmake_minimum_required(VERSION 3.0.0)

project(MatProgramming)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

add_subdirectory(Utils)

add_subdirectory(Function)
add_subdirectory(Point)
add_subdirectory(Matrix)

add_subdirectory(Abstract-solver)
add_subdirectory(Matrix-abstract-solver)

add_subdirectory(Matrix-determinant-solver)
add_subdirectory(Matrix-inverse-solver)
add_subdirectory(Matrix-sylvester-test)

add_subdirectory(Simple-polynomial)
add_subdirectory(Simple-polynomial-real-roots-boundaries-solver)

add_subdirectory(NR)

set(INCLUDE_DIR include)
set(SOURCE_DIR src)

add_library(MatProgramming STATIC ${INCLUDE_DIR}/Header.h ${SOURCE_DIR}/main.cpp)

target_link_libraries(MatProgramming Utils)

target_link_libraries(MatProgramming Point)
target_link_libraries(MatProgramming Matrix)
target_link_libraries(MatProgramming Function)

target_link_libraries(MatProgramming AbstractSolver)
target_link_libraries(MatProgramming MatrixAbstractSolver)
target_link_libraries(MatProgramming MatrixInverseSolver)

target_link_libraries(MatProgramming MatrixDeterminantSolver)

target_link_libraries(MatProgramming MatrixSylvesterTest)

target_link_libraries(MatProgramming SimplePolynomial)
target_link_libraries(MatProgramming SimplePolynomialRealRootsBoundariesSolver)

target_link_libraries(MatProgramming NR)

target_include_directories(MatProgramming PRIVATE ${SOURCE_DIR})
target_include_directories(MatProgramming PUBLIC ${INCLUDE_DIR})