message("BUILD_EXAMPLES is active")
#!------------------------------------------------------
#! ForwardDynamicsCRTP example
add_executable(example_FwdDyn example_FwdDyn.cpp)
target_link_libraries (example_FwdDyn ${GeoMBD_LIBRARIES})
set_target_properties(example_FwdDyn PROPERTIES LINKER_LANGUAGE CXX)
#target_compile_options(example_FwdDyn PRIVATE ${_CXX_FLAGS_RELEASE})
##set_target_properties(example_FwdDyn PROPERTIES COMPILE_FLAGS "-save-temps")
#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example
add_executable(example_FwdDynDiff example_FwdDynDiff.cpp)
target_link_libraries (example_FwdDynDiff ${GeoMBD_LIBRARIES})
set_target_properties(example_FwdDynDiff PROPERTIES LINKER_LANGUAGE CXX)
