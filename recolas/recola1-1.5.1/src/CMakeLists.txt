# File: CMakeLists.txt
# Author: Jean-Nicolas Lang
# Description: CMake configure file for Recola sources
# Last Modified: September 30, 2024

option (static "shared_libraries" OFF)
option (lib_suffix "Whether to add suffix to library" "")
option (rambo "with rambo" OFF)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set (unoptimized "-O0")
if (Fortran_COMPILER_NAME MATCHES "gfortran.*" OR Fortran_COMPILER_NAME MATCHES "f95")
  # gfortran
  set (debug " -g -ffixed-line-length-132")
  set (checks " -fcheck=all -ffpe-trap=invalid,zero,overflow")
  set (inits " -finit-integer=-1000000 -finit-logical=true -finit-real=nan")
  set (warnings " -Wall ")
  set (maths " -fno-fast-math -fno-omit-frame-pointer -ftrapv")
  set (debug "${debug} ${checks} ${inits} ${maths} ${warnings}")
  set (CMAKE_Fortran_FLAGS_DEBUG "${debug}")
  set (CMAKE_Fortran_FLAGS_RELEASE "-ffixed-line-length-132 -fPIC")
  set (optimized "-O2")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (basic " -132")
  set (opt_DEBUG " -g -traceback -O0 -check all")
  set (opt_DEBUG "${opt_DEBUG} -warn all")
  set (CMAKE_Fortran_FLAGS_RELEASE "${basic} ${opt_RELEASE}")
  set (CMAKE_Fortran_FLAGS_DEBUG "${basic} ${opt_DEBUG}")
  set (optimized "-O2 -ipo")
elseif (Fortran_COMPILER_NAME MATCHES "pgf77")
  # g77 (untested)
  set (basic "${basic} -ffixed-line-length-132  -g77libs -Msecond_underscore")
  set (opt_DEBUG "-O0 -Mbounds -frange-check -g -fcheck=all")
  set (opt_DEBUG "${opt_DEBUG} -Wall -Wtabs -Wextra -Wno-unused")
  set (opt_DEBUG "${opt_DEBUG} -Wno-unused-dummy-argument")
  set (opt_DEBUG "${opt_DEBUG} -Wno-unused-parameter")
  set (CMAKE_Fortran_FLAGS_RELEASE "${basic} ${opt_RELEASE}")
  set (CMAKE_Fortran_FLAGS_DEBUG "${basic} ${opt_DEBUG}")
  set (optimized "-O2")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
  set (optimized "-O2")
endif ()

# need preprocessor flag support
set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -cpp")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp")


set(RECOLA_FILES "")

SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_SOURCE_PATH}/recola.f90)

set(RECOLA_CORE_PATH ${RECOLA_SOURCE_PATH}/internal)
set(RECOLA_API_PATH ${RECOLA_SOURCE_PATH})

#========#
#  Core  #
#========#

SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/amplitude_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/collier_interface_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/currents_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/draw_current_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/globals_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/loop_functions_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/loop_vertices_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/model_vertices_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/outgoing_momenta_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/skeleton_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/statistics_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/tables_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/tree_vertices_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/wave_functions_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/wrapper_rcl.f90)

#========#
#  Util  #
#========#

find_path(rambo_path
          NAMES rambo.f90
          PATHS ${RECOLA_CORE_PATH})
find_path(ranlux_path
          NAMES ranlux.F
          PATHS ${RECOLA_CORE_PATH})
if(rambo_path AND ranlux_path)
  message(STATUS "Including rambo in the build.")
  set(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/rambo.f90)
  set(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/ranlux.F)
  add_definitions(-DWITHRAMBO)
else()
  remove_definitions(-DWITHRAMBO)
endif()

IF (with_smtests)
  SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_CORE_PATH}/check_rcl.f90)
ENDIF()


#=======#
#  API  #
#=======#

SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_API_PATH}/input_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_API_PATH}/process_computation_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_API_PATH}/process_definition_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_API_PATH}/process_generation_rcl.f90)
SET(RECOLA_FILES ${RECOLA_FILES} ${RECOLA_API_PATH}/reset_rcl.f90)

string(TOLOWER "${CMAKE_BUILD_TYPE}" CBT)
if (CBT MATCHES debug)
  set(OFLAG ${unoptimized})
else()
  set(OFLAG ${optimized})
endif ()
set_source_files_properties(${RECOLA_FILES} PROPERTIES COMPILE_FLAGS ${OFLAG})


# gather all source files
set(SOURCE ${RECOLA_FILES})

# define the recola library
if (static)
  message(STATUS "Building recola as static library")
  if (lib_suffix)
    set(RECOLA_LIBRARY "recola_${libtype_suffix}")
  else()
    set(RECOLA_LIBRARY "recola")
  endif()
  add_library (${RECOLA_LIBRARY} STATIC ${SOURCE})
else()
  message(STATUS "Building recola as shared library")
  if (lib_suffix)
    set(RECOLA_LIBRARY "recola_${libtype_suffix}")
  else()
    set(RECOLA_LIBRARY "recola")
  endif()
  add_library (${RECOLA_LIBRARY} SHARED ${SOURCE})
endif()
set(RECOLA_LIBRARY ${RECOLA_LIBRARY} PARENT_SCOPE)

#-------------------#
#  Include Collier  #
#-------------------#

SET(COLLIER_MIN_VERSION 1.2.8)

function(reset_collier_cache)
  unset(COLLIER_INCLUDE_DIR CACHE)
  unset(COLLIER_LIBRARY_DIR CACHE)
  unset(COLLIER_LIBRARY_PATH CACHE)
  unset(COLLIER_SYSCONFIG_DIR CACHE)
  unset(collier_DIR CACHE)
endfunction()

# Option 1: collier_path given
if (collier_path)
  unset(COLLIER_SEARCH_PATH CACHE)
  set (COLLIER_SEARCH_PATH ${collier_path} CACHE PATH "Cached path to collier library")
  unset(collier_path CACHE)
  reset_collier_cache()
# Option 2: via environment variable COLLIER_PATH
elseif(DEFINED ENV{COLLIER_PATH})
  unset(COLLIER_SEARCH_PATH CACHE)
  set (COLLIER_SEARCH_PATH $ENV{COLLIER_PATH} CACHE PATH "Cached path to collier library")
  reset_collier_cache()
endif()

find_package(collier ${COLLIER_MIN_VERSION} REQUIRED HINTS ${COLLIER_SEARCH_PATH} NO_DEFAULT_PATH)
set(COLLIER_INCLUDE_DIR ${COLLIER_INCLUDE_DIR} PARENT_SCOPE)
message(STATUS "Collier library path: ${COLLIER_LIBRARY_PATH}")
message(STATUS "Collier include path: ${COLLIER_INCLUDE_DIR}")
message(STATUS "Collier version exact: ${collier_VERSION}")
message(STATUS "RECOLA LIBRARY : ${RECOLA_LIBRARY}")
target_link_libraries (${RECOLA_LIBRARY} ${COLLIER_LIBRARY_PATH})

#----------------------------#
#  Recola file installation  #
#----------------------------#

set(RECOLA_HEADERS "")
set(RECOLA_HEADERS ${RECOLA_HEADERS} "${CMAKE_Fortran_MODULE_DIRECTORY}/recola.mod")
set(RECOLA_HEADERS ${RECOLA_HEADERS} "${RECOLA_INCLUDE_PATH}/recola.h" "${RECOLA_INCLUDE_PATH}/recola_version.h")
set_target_properties(${RECOLA_LIBRARY} PROPERTIES PUBLIC_HEADER "${RECOLA_HEADERS}")

message("")
if (CBT MATCHES debug)
  message("Build type: Debug")
else()
  message("Build type: Release")
endif ()

message ("Recola:  Fortran compiler       = ${Fortran_COMPILER_NAME}")
message ("         Fortran compiler ID    = ${CMAKE_Fortran_COMPILER_ID}")
if (CBT MATCHES "debug")
  message ("         Fortran flags          = ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  message ("         Fortran flags          = ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message ("         Optimization flag      = ${OFLAG}")
message ("         Library install prefix = ${LIB_INSTALL_DIR}")
message ("         Include install prefix = ${INCLUDE_INSTALL_DIR}\n")

include(PackageConfigInstall)
if (static)
  PACKAGE_CONFIG_INSTALL(${RECOLA_LIBRARY} ${RECOLA_VERSION} STATIC)
else()
  PACKAGE_CONFIG_INSTALL(${RECOLA_LIBRARY} ${RECOLA_VERSION} SHARED)
endif()
