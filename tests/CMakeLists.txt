cmake_minimum_required(VERSION 3.19.3)

find_package(Boost COMPONENTS unit_test_framework REQUIRED NO_MODULE)

foreach(test_src IN ITEMS test_lattice test_extended_toric_code_qmc)
    add_executable(${test_src} ${test_src}.cpp)

    target_include_directories(${test_src}
      PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
    )

    target_compile_definitions(${test_src} PRIVATE BOOST_TEST_DYN_LINK)

    target_link_libraries(${test_src}
      PRIVATE
        paratoric::core
        Boost::unit_test_framework
    )

    add_test(NAME ${test_src} COMMAND ${test_src})
endforeach()