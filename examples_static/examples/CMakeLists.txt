message("BUILD_EXAMPLES is active")
#!------------------------------------------------------
#! ForwardDynamicsCRTP examples

add_executable(example_FD_lwr example_FD_lwr.cpp)
target_link_libraries (example_FD_lwr ${GeoMBD_LIBRARIES})
set_target_properties(example_FD_lwr PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_FD_nao example_FD_nao.cpp)
target_link_libraries (example_FD_nao ${GeoMBD_LIBRARIES})
set_target_properties(example_FD_nao PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_FD_hrp2 example_FD_hrp2.cpp)
target_link_libraries (example_FD_hrp2 ${GeoMBD_LIBRARIES})
set_target_properties(example_FD_hrp2 PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_FD_atlas example_FD_atlas.cpp)
target_link_libraries (example_FD_atlas ${GeoMBD_LIBRARIES})
set_target_properties(example_FD_atlas PROPERTIES LINKER_LANGUAGE CXX)

#!------------------------------------------------------
#! ForwardDynamicsDiffCRTP example

add_executable(example_D_FD_lwr example_D_FD_lwr.cpp)
target_link_libraries (example_D_FD_lwr ${GeoMBD_LIBRARIES})
set_target_properties(example_D_FD_lwr PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_D_FD_nao example_D_FD_nao.cpp)
target_link_libraries (example_D_FD_nao ${GeoMBD_LIBRARIES})
set_target_properties(example_D_FD_nao PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_D_FD_hrp2 example_D_FD_hrp2.cpp)
target_link_libraries (example_D_FD_hrp2 ${GeoMBD_LIBRARIES})
set_target_properties(example_D_FD_hrp2 PROPERTIES LINKER_LANGUAGE CXX)

add_executable(example_D_FD_atlas example_D_FD_atlas.cpp)
target_link_libraries (example_D_FD_atlas ${GeoMBD_LIBRARIES})
set_target_properties(example_D_FD_atlas PROPERTIES LINKER_LANGUAGE CXX)

