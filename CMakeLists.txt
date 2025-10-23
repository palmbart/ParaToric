cmake_minimum_required(VERSION 3.23) # for FILE_SET
project(paratoric
  VERSION 1.0
  LANGUAGES CXX C
)

set(PARATORIC_PRE_RELEASE "beta")  # or empty

# Compose full version
set(PARATORIC_VERSION_FULL "${PROJECT_VERSION}")
if(PARATORIC_PRE_RELEASE)  # non-empty
  string(APPEND PARATORIC_VERSION_FULL "-${PARATORIC_PRE_RELEASE}")
endif()

# Use it for packaging/artifacts
set(CPACK_PACKAGE_VERSION "${PARATORIC_VERSION_FULL}")
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")

# Must set CPACK_* before include(CPack)
include(CPack)

# ---- Options ---------------------------------------------------------------

option(PARATORIC_EXPORT_COMPILE_COMMANDS "Export compile_commands.json for tooling" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ${PARATORIC_EXPORT_COMPILE_COMMANDS})

option(PARATORIC_BUILD_PYBIND "Build Python bindings (_paratoric)" OFF)
option(PARATORIC_BUILD_CLI    "Build the paratoric CLI" ${PROJECT_IS_TOP_LEVEL})
option(PARATORIC_BUILD_TESTS  "Build tests" ${PROJECT_IS_TOP_LEVEL})
option(PARATORIC_ENABLE_NATIVE_OPT "Enable CPU-specific optimizations" OFF)
option(PARATORIC_ENABLE_AVX2 "Enable AVX2 flag" OFF)
option(PARATORIC_ENABLE_FAST_MATH "Enable fast-math flag" ON)
option(PARATORIC_LINK_MPI "Try to link to MPI (required on some computing clusters)" OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(paratoric_core STATIC)

# Prefer feature requirements; consumers inherit standard via target
target_compile_features(paratoric_core PUBLIC cxx_std_23)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  if(PARATORIC_ENABLE_NATIVE_OPT)
    target_compile_options(paratoric_core PRIVATE -march=native)
  endif()
  if(PARATORIC_ENABLE_AVX2)
    target_compile_options(paratoric_core PRIVATE -mavx2)
  endif()
  if(PARATORIC_ENABLE_FAST_MATH)
    target_compile_options(paratoric_core PRIVATE -ffast-math)
  endif()
  if(PROJECT_IS_TOP_LEVEL)
    target_compile_options(paratoric_core PRIVATE
      -Wall -Wextra -Wfatal-errors
      -Wno-unused-parameter -Wno-reorder -Wno-unknown-pragmas
      -fdiagnostics-color=always -pipe
    )
  endif()
endif()

set(SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lattice/lattice.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mcmc/extended_toric_code.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mcmc/extended_toric_code_c.cpp
)

set(HDR_PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/paratoric/mcmc/extended_toric_code.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/paratoric/mcmc/extended_toric_code_c.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/paratoric/types/types.hpp
)

# Private headers stay out of install (IDE-visible if you like)
set(HDR_PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lattice/lattice.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mcmc/extended_toric_code_qmc.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rng/rng.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/statistics/autocorrelation.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/statistics/bootstrap.hpp
)

target_sources(paratoric_core
  PRIVATE
    ${SRC_FILES}
  PRIVATE
    FILE_SET private_headers
      TYPE HEADERS
      BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
      FILES ${HDR_PRIVATE}
  PUBLIC
    FILE_SET public_headers
      TYPE HEADERS
      BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
      FILES ${HDR_PUBLIC}
)

target_include_directories(paratoric_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

if (PARATORIC_LINK_MPI)
  find_package(MPI QUIET)
  if(MPI_FOUND)
    target_link_libraries(paratoric_core PRIVATE MPI::MPI_C MPI::MPI_CXX)
  else()
    message(STATUS "MPI not found: building core without MPI")
  endif()
endif()

find_package(Boost COMPONENTS log_setup log REQUIRED NO_MODULE)
target_link_libraries(paratoric_core PUBLIC Boost::log_setup Boost::log)

set_target_properties(paratoric_core PROPERTIES
  OUTPUT_NAME paratoric_core
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME core
)

add_library(paratoric::core ALIAS paratoric_core)

# ---- CLI (optional) --------------------------------------------------------

if(PARATORIC_BUILD_CLI)
  find_package(Boost REQUIRED COMPONENTS program_options log_setup log NO_MODULE)
  find_package(HDF5 REQUIRED COMPONENTS C CXX)

  add_executable(paratoric src/cli/paratoric.cpp)
  target_sources(paratoric
    PRIVATE
      src/io/io.cpp
      src/io/io.hpp
  )
  target_link_libraries(paratoric PRIVATE paratoric::core Boost::program_options Boost::log_setup Boost::log HDF5::HDF5)
  target_include_directories(paratoric
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${CMAKE_CURRENT_SOURCE_DIR}/src
  )
  set_target_properties(paratoric PROPERTIES OUTPUT_NAME paratoric)
  install(TARGETS paratoric RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
  # sometimes useful, symbolic links/copies don't work on macOS due to permission issues
  if(PROJECT_IS_TOP_LEVEL)
    install(TARGETS paratoric RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
  endif()
endif()

# ---- Install & Export ------------------------------------------------------

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS paratoric_core
  EXPORT paratoricTargets
  FILE_SET public_headers
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package name (spell this once and reuse)
set(PARATORIC_PACKAGE_NAME paratoric)

# Generate Config + Version (install-tree)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PARATORIC_PACKAGE_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Ensure the config template exists where we expect it
set(PARATORIC_CONFIG_IN "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PARATORIC_PACKAGE_NAME}Config.cmake.in")
if(NOT EXISTS "${PARATORIC_CONFIG_IN}")
  message(FATAL_ERROR "Missing package config template: ${PARATORIC_CONFIG_IN}")
endif()

configure_package_config_file(
  "${PARATORIC_CONFIG_IN}"
  "${CMAKE_CURRENT_BINARY_DIR}/${PARATORIC_PACKAGE_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PARATORIC_PACKAGE_NAME}"
)

# Install the config + version (install-tree)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PARATORIC_PACKAGE_NAME}Config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PARATORIC_PACKAGE_NAME}ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PARATORIC_PACKAGE_NAME}"
)

# Install the exported targets (install-tree)
install(EXPORT paratoricTargets
  NAMESPACE paratoric::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PARATORIC_PACKAGE_NAME}"
)

# ---- Build-tree package (developer convenience) ---------------------------

# Put build-tree package files in a predictable folder
set(PARATORIC_BUILD_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/cmake")
file(MAKE_DIRECTORY "${PARATORIC_BUILD_PACKAGE_DIR}")

# Export targets for the build tree
export(
  EXPORT paratoricTargets
  NAMESPACE paratoric::
  FILE "${PARATORIC_BUILD_PACKAGE_DIR}/${PARATORIC_PACKAGE_NAME}Targets.cmake"
)

# Generate build-tree Config + Version next to Targets
write_basic_package_version_file(
  "${PARATORIC_BUILD_PACKAGE_DIR}/${PARATORIC_PACKAGE_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${PARATORIC_CONFIG_IN}"
  "${PARATORIC_BUILD_PACKAGE_DIR}/${PARATORIC_PACKAGE_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PARATORIC_PACKAGE_NAME}"
)

# Make development 'find_package(paratoric CONFIG)' work without PATH tweaks
if(PROJECT_IS_TOP_LEVEL AND CMAKE_EXPORT_PACKAGE_REGISTRY)
  export(PACKAGE ${PARATORIC_PACKAGE_NAME})
endif()

# ---- Tests (only when desired) -------------------------------------------

if(PARATORIC_BUILD_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()

# ---- Python bindings (opt-in; subproject-friendly) -----------------------

if(PARATORIC_BUILD_PYBIND)
  add_subdirectory(external/pybind11)
  add_subdirectory(python/bindings)
endif()
