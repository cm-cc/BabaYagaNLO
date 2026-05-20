################################################################################
#                             Recola Configuration                             #
################################################################################

# File: CMakeLists.txt
# Author: Jean-Nicolas Lang
# Description: CMake configure file for Recola
# Last Modified: September 30, 2024

set(CMAKE_MACOSX_RPATH 1) 
set(CMAKE_CACHEFILE_DIR  CMakeFiles/)
cmake_minimum_required (VERSION 3.5...4.0)

#-----------------------#
#  Project declaration  #
#-----------------------#

project (recola Fortran C CXX)

set(RECOLA_VERSION 1.5.1)

#sm testsuite. Checks against Pole & OpenLoops
option(with_smtests "Compile SM testsuite." On)
option(with_all_smtests "Compile the complete SM testsuite." OFF)

set(RECOLA_ROOT_PATH ${PROJECT_SOURCE_DIR})
set(RECOLA_SOURCE_PATH ${RECOLA_ROOT_PATH}/src)
set(RECOLA_INCLUDE_PATH ${RECOLA_ROOT_PATH}/include)

# Set modules directory
set(LIB_LOCAL_DIR ${PROJECT_SOURCE_DIR})
set(SYSCONFIG_LOCAL_DIR ${LIB_LOCAL_DIR}/config)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_LOCAL_DIR}/include)
set(INCLUDE_LOCAL_DIR ${LIB_LOCAL_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_LOCAL_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_LOCAL_DIR})
include_directories (${CMAKE_Fortran_MODULE_DIRECTORY})

# define canonical relative paths
set(LIB_INSTALL_DIR lib CACHE PATH "Relative installation path for the library")
set(INCLUDE_INSTALL_DIR include CACHE PATH "Relative installation path for the headers")
set(SYSCONFIG_INSTALL_DIR lib/cmake CACHE PATH "Relative installation path for the configs")

# Make relative paths absolute (needed later on)
if(NOT IS_ABSOLUTE "${LIB_INSTALL_DIR}")
  set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
endif()
if(NOT IS_ABSOLUTE "${INCLUDE_INSTALL_DIR}")
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
endif()
if(NOT IS_ABSOLUTE "${SYSCONFIG_INSTALL_DIR}")
  set(SYSCONFIG_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${SYSCONFIG_INSTALL_DIR}")
endif()

# set install rpath to lib install dir
SET(CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")

# include check flags (CMake 3.3x feature)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/config)
include(CheckFortranCompilerFlag)

# add recola sources
add_subdirectory(src)

#-----------------------------#
#  Include recola demo files  #
#-----------------------------#

# Define demos directory
set (demos_dir "${RECOLA_ROOT_PATH}/demos")
# Define demo source files
file (GLOB_RECURSE demo_source_files ${demos_dir}/demo*_rcl.f90)
file (GLOB_RECURSE cdemo_source_files ${demos_dir}/cdemo*_rcl.cpp)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*" OR Fortran_COMPILER_NAME MATCHES "f95")
  add_definitions(-DGNU_NAME_MANGLING)
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  add_definitions(-DINTEL_NAME_MANGLING)
endif()

# Create demo executable
foreach (src ${demo_source_files})
  get_filename_component(src_name ${src} NAME_WE)
  add_executable (${src_name} ${demos_dir}/${src_name}.f90)
  set_target_properties ( ${src_name} PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${demos_dir}
                          EXCLUDE_FROM_DEFAULT_BUILD TRUE
                          EXCLUDE_FROM_ALL TRUE           )
  target_link_libraries (${src_name} ${RECOLA_LIBRARY})
endforeach()
foreach (src ${cdemo_source_files})
  get_filename_component(src_name ${src} NAME_WE)
  add_executable (${src_name} ${demos_dir}/${src_name}.cpp)
  set_target_properties ( ${src_name} PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${demos_dir}
                          EXCLUDE_FROM_DEFAULT_BUILD TRUE
                          EXCLUDE_FROM_ALL TRUE           )
  target_link_libraries (${src_name} ${RECOLA_LIBRARY})
endforeach()

#---------#
#  Tests  #
#---------#

if (with_smtests OR with_all_smtests)
  enable_testing()
  add_subdirectory(tests)
endif()
