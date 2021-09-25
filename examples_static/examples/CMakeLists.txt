message("BUILD_EXAMPLES is active")
#!------------------------------------------------------
#! ForwardDynamicsCRTP examples

add_executable(example_FD example_FD.cpp)
target_link_libraries (example_FD ${GeoMBD_LIBRARIES})
set_target_properties(example_FD PROPERTIES LINKER_LANGUAGE CXX)

#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example

add_executable(example_D_FD example_D_FD.cpp)
target_link_libraries (example_D_FD ${GeoMBD_LIBRARIES})
set_target_properties(example_D_FD PROPERTIES LINKER_LANGUAGE CXX)

