cmake_minimum_required(VERSION 3.19.3)

find_package(Boost COMPONENTS unit_test_framework REQUIRED NO_MODULE)

set(test_sources
    lattice/test_lattice.cpp
    mcmc/test_extended_toric_code_qmc.cpp
)

foreach(test_src IN LISTS test_sources)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})

    target_include_directories(${test_name}
      PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
    )

    target_compile_definitions(${test_name} PRIVATE BOOST_TEST_DYN_LINK)

    target_link_libraries(${test_name}
      PRIVATE
        paratoric::core
        Boost::unit_test_framework
    )

    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
