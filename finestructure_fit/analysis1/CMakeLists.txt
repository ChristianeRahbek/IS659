cmake_minimum_required(VERSION 3.16)
project(analysis1)

#set(CMAKE_CXX_STANDARD 14)

#add_executable(analysis1 Analysis.cpp)

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

add_executable(analysis1 Analysis.cpp)
target_link_libraries(analysis1 PRIVATE AUSA::AUSA config++)
target_include_directories(analysis1 PRIVATE ${CMAKE_SOURCE_DIR}/include)
