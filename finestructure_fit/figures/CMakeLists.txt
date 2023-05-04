cmake_minimum_required(VERSION 3.16)
project(figures)

find_package(AUSALIB REQUIRED)
find_package(libconfig++ QUIET)
if(NOT libconfig++_FOUND)
    include(../analysis1/cmake/tryharder.cmake)
endif()

# projectutil.h utility header
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
set(HOME_DIR $ENV{HOME})
configure_file(${CMAKE_SOURCE_DIR}/../analysis1/cmake/projectutil.h.in ${CMAKE_SOURCE_DIR}/../analysis1/include/projectutil.h @ONLY)

link_libraries(AUSA::AUSA config++)
include_directories(include)
include(FetchContent)
FetchContent_Declare(eutil GIT_REPOSITORY https://gitlab.au.dk/ausa/erik/eutil.git)
FetchContent_MakeAvailable(eutil)


add_executable(calibrationProcess calibrationProcess.cpp)
add_executable(Gd148PeakCalibration Gd148PeakCalibration.cpp)
add_executable(SolidAngleCalcAlphaT SolidAngleCalc_alpha_t.cpp)