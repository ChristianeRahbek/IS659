cmake_minimum_required(VERSION 3.15)
project(analysis)

link_libraries(AUSA::AUSA config++)
include_directories(include)

add_executable(bat bat.cpp)
add_executable(bat_split bat_split.cpp)