# cmakelists to compile the quench code
# author: Y.Yang (Kyushu University)
# date: 1st Aug 2016
# modified: 22nd Aug 2016

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# setup project info
project(QuenchCode CXX)
set(QUENCH_MAJOR_VERSION 1)
set(QUENCH_MINOR_VERSION 0)
set(QUENCH_VERSION ${QUENCH_MAJOR_VERSION}.${QUENCH_MINOR_VERSION})
message("-- Quench code version: ${QUENCH_VERSION}")

set(TARGET SimQuench.exe)

# debug
set(CMAKE_VERBOSE_MAKEFILE OFF)

# set c++11 library
set(CMAKE_CXX_COMPILER g++)
add_compile_options(-std=c++11)

# set the path where the executable file should be put
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/bin)

# set the path where the library should be put
set(LIBRARY_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/bin)

# check python
if ( NOT ${PYTHONLIBS_FOUND} )
  message(FATAL_ERROR "python library is not found.")
else ()
  message("-- python version: ${PYTHONLIBS_VERSION_STRING}")
  message("-- python library: ${PYTHON_LIBRARIES}")
endif ( NOT ${PYTHONLIBS_FOUND} )

# find FindROOT.cmake
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)
if ( NOT ROOT_FOUND )
  message(FATAL_ERROR "package is required ROOT, but not found it.")
endif ( NOT ROOT_FOUND )

# check the version of ROOT, currently the required version is 5.34
if ( ${ROOT_VERSION} MATCHES "^5\\.34" )
  message("-- ROOT version: ${ROOT_VERSION}")
else ()
  message("-- ROOT version: ${ROOT_VERSION}")
  message(FATAL_ERROR "quench code is not supporting this version of ROOT.") 
endif ()

# append the root library
set(libs)
list(APPEND libs ${ROOT_LIBRARIES})
include_directories(${ROOT_INCLUDE_DIR})

# include quench inculde files
include_directories(include)
add_subdirectory(src)

add_executable(${TARGET} main.cpp)
target_link_libraries(${TARGET} ${libs} SimQuench)
