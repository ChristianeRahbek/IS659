cmake_minimum_required(VERSION 3.15)
project(is659)

# project dependencies
find_package(AUSALIB REQUIRED) # provides ROOT
find_package(libconfig++)
#if(NOT libconfig++_FOUND)
#  include(cmake/pkgconfig_libconfig.cmake)
#else()
#  message(FATAL_ERROR "CMake config file for libconfig++-dev found.")
#endif()

# fetch commonly used project utilities for convenience
file(READ ${CMAKE_SOURCE_DIR}/VERSION CMAKE_PROJECT_VERSION)
include(FetchContent)
FetchContent_Declare(cmake_common GIT_REPOSITORY https://gitlab.au.dk/ausa/cmake_common.git)
FetchContent_MakeAvailable(cmake_common)
FetchContent_Declare(eutil GIT_REPOSITORY https://gitlab.au.dk/ausa/erik/eutil.git)
FetchContent_MakeAvailable(eutil)
include_directories(${CMAKE_SOURCE_DIR}/include)

# locate unpacked data dir and make it available
include(cmake/path_config.cmake)

# subprojects defined in subdirectories' CMakeLists.txt files
include(cmake/add_subproject.cmake)
list(APPEND subprojects analysis calibration sorting)
add_subproject("${subprojects}")