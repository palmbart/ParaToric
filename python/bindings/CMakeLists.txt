# /python/bindings/CMakeLists.txt

# Assumes:
# - Top-level already did: add_subdirectory(external/pybind11)
# - Target paratoric_core exists
# - C++ standard is set globally (C++23 in your top-level)

find_package(Python3 COMPONENTS Development.Module NumPy REQUIRED)

set(BINDINGS_PY_PACKAGE_DIR
    "${PROJECT_SOURCE_DIR}/python/paratoric"
    CACHE PATH "Directory to place the built Python extension")

set(BINDINGS_MODULE_NAME "_paratoric" CACHE STRING "Python extension module name")

pybind11_add_module(${BINDINGS_MODULE_NAME}
    py_extended_toric_code.cpp
)

target_include_directories(${BINDINGS_MODULE_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(${BINDINGS_MODULE_NAME}
    PRIVATE
        paratoric::core
        pybind11::module
        Python3::NumPy
)

set_property(TARGET paratoric_core PROPERTY POSITION_INDEPENDENT_CODE ON)

set_target_properties(${BINDINGS_MODULE_NAME} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${BINDINGS_PY_PACKAGE_DIR}"
    # For MSVC/Windows:
    RUNTIME_OUTPUT_DIRECTORY "${BINDINGS_PY_PACKAGE_DIR}"
)

target_compile_definitions(${BINDINGS_MODULE_NAME} PRIVATE
    $<$<CONFIG:Debug>:PYBIND11_DETAILED_ERROR_MESSAGES=1>
)

message(STATUS "Python extension '${BINDINGS_MODULE_NAME}' will be written to: ${BINDINGS_PY_PACKAGE_DIR}")

set(_PY_SITELIB "")
if (DEFINED Python3_SITELIB AND Python3_SITELIB)
  set(_PY_SITELIB "${Python3_SITELIB}")
elseif (DEFINED Python3_SITEARCH AND Python3_SITEARCH)
  set(_PY_SITELIB "${Python3_SITEARCH}")
else()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import site,sys; print(next((p for p in site.getsitepackages() if 'site-packages' in p), site.getusersitepackages()))"
    OUTPUT_VARIABLE _PY_SITELIB
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

# 1) Optional: copy into current Python's site-packages
if (PARATORIC_INSTALL_TO_SITE)
  message(STATUS "On install, will copy '${BINDINGS_MODULE_NAME}' to: ${_PY_SITELIB}/paratoric")

  install(TARGETS ${BINDINGS_MODULE_NAME}
          LIBRARY DESTINATION "${_PY_SITELIB}/paratoric"
          RUNTIME DESTINATION "${_PY_SITELIB}/paratoric")

  install(FILES
      "${PROJECT_SOURCE_DIR}/python/paratoric/__init__.py"
      "${PROJECT_SOURCE_DIR}/python/paratoric/_paratoric.pyi"
      "${PROJECT_SOURCE_DIR}/python/paratoric/py.typed"
      DESTINATION "${_PY_SITELIB}/paratoric")
endif()

# 2) Optional: run pip editable install (package metadata/entry-points)
if (PARATORIC_PIP_EDITABLE_INSTALL)
  if (PARATORIC_PIP_OFFLINE_INSTALL)
    install(CODE
      "execute_process(COMMAND \"${Python3_EXECUTABLE}\" -m pip install -e \"${PROJECT_SOURCE_DIR}/python\" --no-deps --no-index --no-build-isolation \
                      RESULT_VARIABLE ec)
      if (NOT ec EQUAL 0)
        message(FATAL_ERROR \"pip editable install failed (exit code: \${ec})\")
      endif()"
    )
  else()
    install(CODE
      "execute_process(COMMAND \"${Python3_EXECUTABLE}\" -m pip install -e \"${PROJECT_SOURCE_DIR}/python\" \
                      RESULT_VARIABLE ec)
      if (NOT ec EQUAL 0)
        message(FATAL_ERROR \"pip editable install failed (exit code: \${ec})\")
      endif()"
    )
  endif()
  
endif()