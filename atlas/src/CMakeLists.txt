# cmakelists to compile the quench code
# author: Y.Yang (Kyushu University)
# date: 1st Aug 2016
# modifed: 22nd Aug 2016

cmake_minimum_required(VERSION 2.8)

set(QUENCH_STATIC_LIB SimAtlasQuench)

file(GLOB QUENCH_SRCS "*.cxx")

add_library(${QUENCH_STATIC_LIB} STATIC ${QUENCH_SRCS})
