cmake_minimum_required(VERSION 3.16)

project(groopsio LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 14)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Python3_NumPy_INCLUDE_DIRS})

find_package(EXPAT  REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

find_library(LIB_Z z)
if(LIB_Z)
  find_path(ZLIB_INCLUDE_DIR NAMES zlib.h)
  include_directories(${ZLIB_INCLUDE_DIR})
  set(BASE_LIBRARIES ${BASE_LIBRARIES} ${LIB_Z})
else()
  add_definitions(-DNOLIB_Z)
  message(WARNING "Z library *NOT* found (https://www.zlib.net). groopsio will not be able to read/write compressed *.gz files.")
endif()

add_library(lapack_lite UNKNOWN IMPORTED)
set_property(TARGET lapack_lite PROPERTY IMPORTED_LOCATION "${LAPACK_LITE}")

set(BASE_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXPAT_LIBRARIES} ${LIB_Z} ${Python3_LIBRARIES} lapack_lite stdc++fs)

include(sources.list)
set(GROOPS_SOURCE_DIR "" CACHE PATH "path to GROOPS source directory")
include_directories(${GROOPS_SOURCE_DIR})
list(TRANSFORM SOURCES PREPEND ${GROOPS_SOURCE_DIR}//)

set(EXTENSION_LIBRARY_NAME "groopsiobase" CACHE STRING "name of the Python extension library")
add_library(${EXTENSION_LIBRARY_NAME} SHARED ${PROJECT_SOURCE_DIR}/src/groopsio.cpp ${SOURCES})

SET_TARGET_PROPERTIES(${EXTENSION_LIBRARY_NAME} PROPERTIES SUFFIX "")
SET_TARGET_PROPERTIES(${EXTENSION_LIBRARY_NAME} PROPERTIES PREFIX "")

target_link_libraries(${EXTENSION_LIBRARY_NAME} ${BASE_LIBRARIES})
