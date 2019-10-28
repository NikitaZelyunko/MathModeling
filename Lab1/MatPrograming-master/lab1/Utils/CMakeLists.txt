cmake_minimum_required(VERSION 3.0.0)

project(Utils)

set(SOURCE_DIR src)
set(INCLUDE_DIR include)

add_library(Utils STATIC ${INCLUDE_DIR}/Utils.h ${SOURCE_DIR}/Utils.cpp)

# Add visible include directories
# Now include "Utils.h" not throw error
target_include_directories(Utils PRIVATE ${SOURCE_DIR})
target_include_directories(Utils PUBLIC ${INCLUDE_DIR})