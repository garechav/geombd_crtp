cmake_minimum_required(VERSION 2.8.3)
project(Examples_Static)

set(PACKAGE_AUTHOR "Alvaro Paz, Gustavo Arechavaleta")
set(PACKAGE_MAINTAINER "Alvaro Paz")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-msse -msse2 -msse3 -mavx -march=native -O3 -pedantic -DNDEBUG")

list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../geombd/cmake/linux-64)
find_package(geombd)
include_directories(${GeoMBD_INCLUDE_DIRS})

set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../eigen3")

include_directories(${EIGEN3_INCLUDE_DIR})
message(${EIGEN3_INCLUDE_DIR})

add_subdirectory(examples)
