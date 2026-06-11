# File: CMakeLists.txt
# Author: Jean-Nicolas Lang
# Description: CMake configure file for tests
# Last Modified: June 26, 2025

set(CMAKE_CTEST_COMMAND ctest) 
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
include_directories(${COLLIER_INCLUDE_DIR})

function(new_fortran_test testname testsource)
  add_executable(${testsource} EXCLUDE_FROM_ALL ${testsource}.f90)
  target_link_libraries (${testsource} recola)
  add_test(${testname} ${testsource})
  add_dependencies(check ${testsource})
endfunction(new_fortran_test)

function(new_cpp_test testname testsource)
    add_executable(${testsource} EXCLUDE_FROM_ALL ${testsource}.cpp)
    target_link_libraries (${testsource} recola)
    add_test(${testname} ${testsource})
    add_dependencies(check ${testsource})
endfunction(new_cpp_test)

if (with_all_smtests)
    new_fortran_test("two_to_two_processes" twoTOtwo_all)
    new_fortran_test("two_to_three_processes" twoTOthree_all)
    new_fortran_test("two_to_four_processes" twoTOfour_all)
else()
    new_fortran_test("two_to_two_processes" twoTOtwo)
    new_fortran_test("two_to_three_processes" twoTOthree)
    new_fortran_test("two_to_four_processes" twoTOfour)
endif()
new_fortran_test("helicity_projection" helicity_projection_soft_limit)
new_fortran_test("momenta_correction" momenta_correction)
new_cpp_test("cpp_interface" cinterface)
new_cpp_test("helicity_projection_pole_approx_vs_onshell" helicity_projection_pole_vs_onshell)
