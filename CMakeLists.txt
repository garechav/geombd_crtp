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


#-mavx -march=native -O3
#avx,avx2,fma
set(CMAKE_CXX_FLAGS_RELEASE "-msse -msse2 -msse3 -mavx -march=native -O3 -pedantic -DNDEBUG")


#set(CMAKE_CXX_FLAGS_RELEASE "-mavx -march=native -O3 -ftree-vectorize -msse -msse2 -msse3 -fopt-info-vec-missed")
#set(CMAKE_CXX_FLAGS_RELEASE "-mavx -march=native -O3 -flto -unroll-aggressive -finline-atomics -ffast-math -freciprocal-math")


include_directories("include")

set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
message(${EIGEN3_INCLUDE_DIR})

#!------------------------------------------------------

add_subdirectory(include/geombd/io)

#!------------------------------------------------------
#! Under development
#!------------------------------------------------------
#! ForwardDynamicsCRTP
add_library(forward-dynamics-crtp SHARED ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/ForwardDynamicsCRTP.hpp)
set_target_properties(forward-dynamics-crtp PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(forward-dynamics-crtp mysincos)
#!------------------------------------------------------
#! ForwardDynamicsCRTP example
add_executable(example_CRTP_05 examples/example_CRTP_05.cpp)
target_link_libraries (example_CRTP_05 forward-dynamics-crtp ${PACKAGE_NAME}-io)# geombd-core geombd-dynamics)
#target_compile_options(example_CRTP_05 PRIVATE ${_CXX_FLAGS_RELEASE})
##set_target_properties(example_CRTP_05 PROPERTIES COMPILE_FLAGS "-save-temps")
#!------------------------------------------------------
#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP
add_library(forward-dynamics-diff-crtp SHARED ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/DForwardDynamicsCRTP.hpp)
#target_compile_options(forward-dynamics-diff-crtp PRIVATE ${_CXX_FLAGS})
set_target_properties(forward-dynamics-diff-crtp PROPERTIES LINKER_LANGUAGE CXX)
#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example
add_executable(example_CRTP_06 examples/example_CRTP_06.cpp)
target_link_libraries (example_CRTP_06 forward-dynamics-diff-crtp ${PACKAGE_NAME}-io)

#!------------------------------------------------------
#! Parser example
#add_executable(example_parser examples/example_parser.cpp)
#target_link_libraries (example_parser geombd-core ${PACKAGE_NAME}-io forward-dynamics-diff-crtp geombd-dynamics)
#!------------------------------------------------------


# Toggle comment for creating the static library of geombd
#--------------------------------------------------------
file(GLOB geoParser ${PACKAGE_NAME}-io)
file(GLOB fwdDyn ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/ForwardDynamicsCRTP.hpp)
file(GLOB fwdDynDiff ${PROJECT_SOURCE_DIR}/include/geombd/CRTP/DForwardDynamicsCRTP.hpp)

ADD_LIBRARY( geombd STATIC ${fwdDyn} ${fwdDynDiff} ${geoParser} )
set_target_properties(geombd PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries( geombd
                      ${EIGEN3_LIBRARIES})
#--------------------------------------------------------
