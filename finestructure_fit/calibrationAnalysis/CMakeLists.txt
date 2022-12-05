cmake_minimum_required(VERSION 3.16)
project(calibrationAnalysis)

#From Erik's file

find_package(AUSALIB REQUIRED)
find_package(libconfig++ QUIET)
if(NOT libconfig++_FOUND)
    include(cmake/tryharder.cmake)
endif()

# projectutil.h utility header
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
set(HOME_DIR $ENV{HOME})
configure_file(${CMAKE_SOURCE_DIR}/cmake/projectutil.h.in ${CMAKE_SOURCE_DIR}/include/projectutil.h @ONLY)

add_executable(calibrationAnalysis Analysis.cpp)
target_link_libraries(calibrationAnalysis PRIVATE AUSA::AUSA config++)
target_include_directories(calibrationAnalysis PRIVATE ${CMAKE_SOURCE_DIR}/include)
