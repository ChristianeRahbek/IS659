# 2021: Debian, and hence Ubuntu, still does not provide a version of libconfig++ with a Cmake config file.
# Other systems do. When Debian gets onboard with the rest of us, this check can probably be safely removed.
message(STATUS "CMake config file for libconfig++-dev not found.")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
message(STATUS "Trying with PkgConfig")
find_package(PkgConfig QUIET)
pkg_search_module(libconfig++ QUIET libconfig++)
if(libconfig++_FOUND)
  message(STATUS "Trying with PkgConfig - success")
else()
  message(STATUS "Trying with PkgConfig - failed")
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)