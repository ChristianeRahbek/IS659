cmake_minimum_required(VERSION 3.16)
project(analysis1)

find_package(AUSALIB REQUIRED)
find_package(libconfig++ QUIET)
if(NOT libconfig++_FOUND)
    include(cmake/tryharder.cmake)
endif()

# projectutil.h utility header
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
set(HOME_DIR $ENV{HOME})
configure_file(${CMAKE_SOURCE_DIR}/cmake/projectutil.h.in ${CMAKE_SOURCE_DIR}/include/projectutil.h @ONLY)

link_libraries(AUSA::AUSA config++)
include_directories(include)
include(FetchContent)
FetchContent_Declare(eutil GIT_REPOSITORY https://gitlab.au.dk/ausa/erik/eutil.git)
FetchContent_MakeAvailable(eutil)


add_executable(analysis1 Analysis.cpp)
#add_executable(bat_split bat_split.cpp)
add_executable(hstack hstack.cpp)
add_executable(AngleAnalysis AngleAnalysis.cpp)
add_executable(BetaAnalysis Analysis_with_betas.cpp)
add_executable(BetaCoincidenceAnalysis BetaCoincidenceAnalysis.cpp)
add_executable(ThreeCoincidenceAnalysis ThreeCoincidenceAnalysis.cpp)
add_executable(TwoCoincidenceAnalysis TwoCoincidenceAnalysis.cpp)
#add_executable(KolmogorovTest KolmogorovTest.cpp)
#target_link_libraries(analysis1 PRIVATE AUSA::AUSA config++)
#target_include_directories(analysis1 PRIVATE ${CMAKE_SOURCE_DIR}/include)
#target_link_libraries(bat_split PRIVATE AUSA::AUSA config++)
#target_include_directories(bat_split PRIVATE ${CMAKE_SOURCE_DIR}/include)
