## CMakeLists to compile the example code
## Author: Y.Yang

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

## setup project info
project(SimAtlasMagQuench CXX)

set(TARGET SimCometMagQuench.exe)

set(QUENCH_ROOT $ENV{SOFTWARE}/cometscqc)

## setup compiler and c++11
set(CMAKE_CXX_COMPILER g++)
add_compile_options(-std=c++11)

## find ROOT
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/root/cmake)
find_package(ROOT REQUIRED)
if ( NOT ROOT_FOUND )
  message(FATAL_ERROR "package is required ROOT, but not found it.")
endif ( NOT ROOT_FOUND )

## append root library
set(libs)
list(APPEND libs ${ROOT_LIBRARIES})
include_directories(${ROOT_INCLUDE_DIR})

## append simquench library
list(APPEND libs ${QUENCH_ROOT}/build/bin/libSimQuench.a)
include_directories(${QUENCH_ROOT}/include)

## append simcometquench library
include_directories(include)
add_subdirectory(src)

add_executable(${TARGET} main.cxx)
target_link_libraries(${TARGET} SimAtlasQuench ${libs})
