#
# Copyright (c) 2021 Cinvestav
#
cmake_minimum_required(VERSION 3.10.2)
#SET(CMAKE_EXPORT_COMPILE_COMMANDS ON) # json compile flags
message("Module geombd-io")
#file(GLOB inc_io ${INC_FOLDER}/io/*.h ${INC_FOLDER}/io.h)
file(GLOB inc_io *.hpp)

set(CMAKE_CXX_STANDARD 17)

add_library(${PACKAGE_NAME}-io SHARED parser.cpp ${inc_io})
target_compile_options(${PACKAGE_NAME}-io PRIVATE ${_CXX_FLAGS})
