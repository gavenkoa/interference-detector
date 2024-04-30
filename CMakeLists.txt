cmake_minimum_required(VERSION 3.5)
project(interference-detector)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# message("Hello world! " $ENV{CONDA_LIB_ROOT} " some")
include(CMakePrintHelpers)
# cmake_print_variables(ENV{CONDA_LIB_ROOT})
# cmake_print_variables(PROJECT_SOURCE_DIR)

# https://cmake.org/cmake/help/latest/command/if.html
if (DEFINED ENV{CONDA_LIB_ROOT})
  if (EXISTS $ENV{CONDA_LIB_ROOT})
    include_directories(${PROJECT_NAME} AFTER $ENV{CONDA_LIB_ROOT}/include)
    link_directories(${PROJECT_NAME} AFTER $ENV{CONDA_LIB_ROOT}/lib)
    message("Added CONDA_LIB_ROOT to INCLUDE / LIB path: " $ENV{CONDA_LIB_ROOT})
  else()
    message("Broken path to CONDA_LIB_ROOT : " $ENV{CONDA_LIB_ROOT})
  endif()
endif()

add_executable(interference-detector src/main.c)
# cmake_print_variables(CMAKE_FIND_LIBRARY_PREFIXES)
# set(CMAKE_STATIC_LIBRARY_PREFIX "")
# set_target_properties(interference-detector PROPERTIES PREFIX "")
target_link_libraries(interference-detector rtlsdr fftw3 m)

