cmake_minimum_required(VERSION 3.0.0)

project(AbstractSolver)

set(SOURCE_DIR src)
set(INCLUDE_DIR include)

add_library(AbstractSolver STATIC 
    ${INCLUDE_DIR}/AbstractSolver.h
    ${SOURCE_DIR}/AbstractSolver.cpp
)

target_include_directories(AbstractSolver PRIVATE ${SOURCE_DIR})
target_include_directories(AbstractSolver PUBLIC ${INCLUDE_DIR})