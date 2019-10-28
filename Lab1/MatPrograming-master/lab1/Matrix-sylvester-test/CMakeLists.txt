cmake_minimum_required(VERSION 3.0.0)

set(THIS_PROJECT_NAME MatrixSylvesterTest)

project(${THIS_PROJECT_NAME})

set(INCLUDE_DIR include)
set(SRC_DIR src)

add_library(${THIS_PROJECT_NAME} STATIC
        ${INCLUDE_DIR}/${THIS_PROJECT_NAME}.h
        ${SRC_DIR}/${THIS_PROJECT_NAME}.cpp)

target_link_libraries(${THIS_PROJECT_NAME} Matrix)
target_link_libraries(${THIS_PROJECT_NAME} MatrixAbstractSolver)
target_link_libraries(${THIS_PROJECT_NAME} Point)

target_include_directories(${THIS_PROJECT_NAME} PRIVATE ${SRC_DIR})
target_include_directories(${THIS_PROJECT_NAME} PUBLIC ${INCLUDE_DIR})