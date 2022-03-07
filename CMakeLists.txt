cmake_minimum_required(VERSION 3.16)

project(groopsio LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 14)

find_package(Python COMPONENTS Development NumPy)
include_directories(${Python_INCLUDE_DIRS})
include_directories(${Python_NumPy_INCLUDE_DIRS})

find_package(EXPAT  REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

find_package(BLAS REQUIRED)

find_library(LIB_Z z)
if(LIB_Z)
  find_path(ZLIB_INCLUDE_DIR NAMES zlib.h)
  include_directories(${ZLIB_INCLUDE_DIR})
else()
  add_definitions(-DNOLIB_Z)
  message(WARNING "Z library *NOT* found (https://www.zlib.net). groopsio will not be able to read/write compressed *.gz files.")
endif()

set(BASE_LIBRARIES ${EXPAT_LIBRARIES} ${BLAS_LIBRARIES} ${LIB_Z} ${Python_LIBRARIES} stdc++fs)

include(sources.list)
set(GROOPS_SOURCE_DIR "" CACHE PATH "path to GROOPS source directory")
include_directories(${GROOPS_SOURCE_DIR})
list(TRANSFORM SOURCES PREPEND ${GROOPS_SOURCE_DIR}//)

set(EXTENSION_LIBRARY_NAME "groopsiobase" CACHE STRING "name of the Python extension library")
add_library(${EXTENSION_LIBRARY_NAME} SHARED ${PROJECT_SOURCE_DIR}/src/groopsio.cpp ${SOURCES})

SET_TARGET_PROPERTIES(${EXTENSION_LIBRARY_NAME} PROPERTIES SUFFIX "")
SET_TARGET_PROPERTIES(${EXTENSION_LIBRARY_NAME} PROPERTIES PREFIX "")

target_link_libraries(${EXTENSION_LIBRARY_NAME} ${BASE_LIBRARIES})
