cmake_minimum_required(VERSION 3.0.0)

project(NR)

set(INCLUDE_DIR include)
set(SOURCE_DIR src)

add_library(NR STATIC 
    ${INCLUDE_DIR}/NR.h
    ${SOURCE_DIR}/NR.cpp
)

target_include_directories(NR PRIVATE ${SOURCE_DIR})
target_include_directories(NR PUBLIC ${INCLUDE_DIR})

target_link_libraries(NR Utils)
