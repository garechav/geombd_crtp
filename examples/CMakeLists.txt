message("BUILD_EXAMPLES is active")

#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example
add_executable(example_FwdDynDiff example_FwdDynDiff.cpp)
target_link_libraries (example_FwdDynDiff forward-dynamics-diff-crtp ${PACKAGE_NAME}-io)
