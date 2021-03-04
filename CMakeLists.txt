cmake_minimum_required(VERSION 3.17)
project(GaussMethod)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(MathMatrix)

add_library(${PROJECT_NAME} GaussMethod.cpp GaussMethod.h)

target_include_directories(${PROJECT_NAME} PUBLIC MathMatrix)

target_link_libraries(${PROJECT_NAME} MathMatrix)