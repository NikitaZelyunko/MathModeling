cmake_minimum_required(VERSION 3.0.0)

project(Point)

set (INCLUDE_DIR include)
set (SOURCE_DIR src)

add_library(Point STATIC ${SOURCE_DIR}/Point.cpp ${INCLUDE_DIR}/Point.h)

target_link_libraries(Point Utils)

target_include_directories(Point PRIVATE ${SOURCE_DIR})
target_include_directories(Point PUBLIC ${INCLUDE_DIR})
