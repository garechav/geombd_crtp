message("BUILD_EXAMPLES is active")

#!------------------------------------------------------
#! ForwardDynamicsCRTP example
add_executable(example_FwdDyn example_FwdDyn.cpp)
target_link_libraries (example_FwdDyn forward-dynamics-crtp ${PACKAGE_NAME}-io)
#target_compile_options(example_FwdDyn PRIVATE ${_CXX_FLAGS_RELEASE})
##set_target_properties(example_FwdDyn PROPERTIES COMPILE_FLAGS "-save-temps")

#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example
add_executable(example_FwdDynDiff example_FwdDynDiff.cpp)
target_link_libraries (example_FwdDynDiff forward-dynamics-diff-crtp ${PACKAGE_NAME}-io)

#!------------------------------------------------------
#! Forward Dynamics + Diff example
add_executable(example_FwdDyn_plus_Diff example_FwdDyn_plus_Diff.cpp)
target_link_libraries (example_FwdDyn_plus_Diff forward-dynamics-crtp forward-dynamics-diff-crtp ${PACKAGE_NAME}-io)

