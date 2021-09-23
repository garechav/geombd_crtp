cmake_minimum_required(VERSION 2.8.3)
project(GeoMBD)

set(PACKAGE_NAME "geombd")
set(PACKAGE_VERSION "1.0.0")
set(PACKAGE_SO_VERSION "1.0")
set(PACKAGE_DESCRIPTION "Geometric Multibody Dynamics Library")
set(PACKAGE_AUTHOR "Alvaro Paz, Gustavo Arechavaleta")
set(PACKAGE_MAINTAINER "Alvaro Paz, Gustavo Arechavaleta")

set(SRC_FOLDER  ${PROJECT_SOURCE_DIR}/src)
set(TEST_FOLDER ${PROJECT_SOURCE_DIR}/test)
set(INC_FOLDER  ${PROJECT_SOURCE_DIR}/include/${PACKAGE_NAME})

#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#SET(CMAKE_EXPORT_COMPILE_COMMANDS ON) # json compile flags

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_EXTENSIONS ON)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)


#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-msse -msse2 -msse3")


set(CMAKE_CXX_FLAGS_RELEASE "-msse -msse2 -msse3 -mavx -march=native -O3 -pedantic -DNDEBUG")


include_directories("include")

set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
message(${EIGEN3_INCLUDE_DIR})

#!------------------------------------------------------

add_subdirectory(include/geombd/io)

#!------------------------------------------------------
#! ForwardDynamicsCRTP
add_library(forward-dynamics-crtp SHARED ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/ForwardDynamicsCRTP.hpp)
set_target_properties(forward-dynamics-crtp PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(forward-dynamics-crtp mysincos)

#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP
add_library(forward-dynamics-diff-crtp SHARED ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/DForwardDynamicsCRTP.hpp)
#target_compile_options(forward-dynamics-diff-crtp PRIVATE ${_CXX_FLAGS})
set_target_properties(forward-dynamics-diff-crtp PROPERTIES LINKER_LANGUAGE CXX)


add_subdirectory(examples)


## Toggle comment for creating the static library of geombd
##--------------------------------------------------------
#file(GLOB src_io ${INC_FOLDER}/io/*.cpp)
##file(GLOB geoParser ${PACKAGE_NAME}-io)
##file(GLOB fwdDyn ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/ForwardDynamicsCRTP.hpp)
##file(GLOB fwdDynDiff ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/DForwardDynamicsCRTP.hpp)

#ADD_LIBRARY( geombd STATIC ${src_io} )
#set_target_properties(geombd PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries( geombd
#                      ${EIGEN3_LIBRARIES})
###--------------------------------------------------------
