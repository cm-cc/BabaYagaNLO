# RECOLA README #

* * *

[**Recola**][1] is a Fortran95 computer program for the automated
generation and numerical computation of EW and QCD amplitudes in the Standard
Model at next-to-leading order.

* * *

## Features

+ Amplitudes in the 't Hooft-Feynman gauge

+ Complex-mass scheme for unstable particles implemented

+ Possibility of consistent isolation of resonant contributions and of applying the pole approximation

+ Dimensional regularization for ultraviolet and infrared singularities, with the possibility of treating collinear and soft singularities in mass regularization

+ Various renormalization schemes supported for the electromagnetic coupling constant

+ Dynamical Nf-flavour renormalization scheme for the strong coupling constant

+ Computation of next-to-leading-order amplitudes for all helicities and colour structures

+ Computation of next-to-leading-order squared amplitudes summed/averaged over spin and colour

+ Computation of colour- and/or spin-correlated leading-order squared amplitudes for dipole subtraction

+ Computation of colour- and/or spin-correlated leading-order squared amplitudes for dipole subtraction

+ Polarizations selection for internal massive fermions and vector bosons

* * *

## Installation

### External dependencies

+ [Collier][2] tensor integral library
+ [CMake][3] build system

### Summary of set up

This *Recola* package is a standalone version and requires the user to resolve dependencies to *Collier* by hand.  We assume the [Collier][2] library has been built as shared libraries and ready to get linked to *Recola*. The compilation of *Recola* consists of the following steps:

1. Configuration of *Recola* with `cmake`
2. Compilation of *Recola* with `make`

### **QUICK** configuration and installation of **Recola**

Obtain a copy of *Recola* from <http://recola.hepforge.org/> and extract it via the shell command

    tar -zxvf recola-X.Y.Z.tar.gz

where X.Y.Z is the current *Recola* version. Then, switch to the build directory and run

    cd recola-X.Y.Z/build
    cmake [options] .. -Dcollier_path=<path_to_collier>
    make [options]

where \<path_to_collier\> points to the directory containing the compiled *Collier* library.

* * *


### Running demo files

The [Recola][1] demo files are located in:

    recola-X.Y.Z/demos

The Fortran and C++ demo files can be run by invoking the run script

    ./run <demofile>

More information on the demo files can be found in the official manual.

* * *

### Summary of CMake and Makefile options

|      CMake option      | Value         | Short description  |
| :--------------------: | :------------ |:------------------ |
| collier path           | Path          | Absolute or relative path to the *Collier* library. |
| static                 | On/Off        | Compile the library as a shared or static library. |
| with_smtests           | On/Off        | Run tests against *Pole* and *OpenLoops*. |
| CMAKE_BUILD_TYPE       | Debug/Release | Set the compiler flags. By default Release flags (optimized) are used. |
| CMAKE_Fortran_COMPILER | Path/Name     | Set the Fortran compiler either via executable name or the absolute path to executable. |
| CMAKE_INSTALL_PREFIX   | Path          | Set the installation prefix. |

|  Makefile option | Value      | Short description  |
| :--------------: | :--------- |:------------------ |
| -j               | Integer    | Number of threads for compilation. |
| VERBOSE          | True/False | Enable *verbose* compilation. In this mode all compilation flags are visible to the user. |

* * *

### Running testsuite

The test suite can be enabled via:

    cmake [options] .. -Dwith_smtests=On -Dcollier_path=<path_to_collier>

After (re-)compilation the testruns can be initiated via

    make check

in the recola build directory.

* * *

## General instruction when using Recola in a Fortran/C++ program

In order to use Recola its modules have to be loaded:

+ Fortran: `use recola`

+ C++: `#include "recola.h"`

in the preamble of the respective code, and the library `librecola.so`(or `librecola.a`, `librecola.dynlib`) has to be supplied to the linker. This gives access to the public functions and subroutines of the Recola library. The names of all these routines end with the suffix `_rcl`.

Typically, an application of Recola involves the following five steps:

   - **Step 1**: *Setting input parameters (optional)*

     The input needed for the computation of processes can be set by calling dedicated subroutines as provided by Recola. See Section 4.2 *Input subroutines* in the official manual for the full set of supported methods.  Since Recola provides default values for all input parameters, this first step is optional.

   - **Step 2**: *Process definition*

     The calculation of matrix elements for one or more processes requires each process to be declared and labelled with a unique identifier. This is done by calling the `define_process_rcl`.

   - **Step 3**: *Process generation*

     In the next step the subroutine `generate_processes_rcl` is called which triggers the initialization of the complete list of processes defined in step 2.

   - **Step 4**: *Process computation*

     Recola is now ready to calculate amplitudes for any of the processes defined in step 2. The computation of the amplitude and of the squared amplitude is performed by means of the subroutine `compute_process_rcl`, which uses the process-dependent information derived in step 3. The subroutine `compute_process_rcl` is called with the momenta of the external particles provided by the user.

   - **Step 5**: *Resetting Recola*

     Finally, by calling the subroutine `reset_recola_rcl`, the process-dependent information generated in steps 2-4 is deleted and the corresponding memory is deallocated. The input variables keep their values defined in step 1 before.

Note that these steps have to be followed in the order given above.  In particular, after step 3 no additional processes can be defined unless Recola is reset (step 5). After step 5 the user can restart with step 1 or step 2.

* * *

## Contributing

We use a branch-based workflow with `testing` and `stable` as integration branches.

### External contributors

External contributors should open a merge request (MR) targeting the `testing` branch.

### Internal contributors and releases

Internal development and releases follow this flow:

1. Work on an issue.
2. Open an MR targeting `testing`.
3. After validation, open an MR from `testing` to `stable`.
4. Create the release from `stable` and publish it with a tag.

* * *

[1]: https://recola.gitlab.io/  "RECOLA"
[2]: http://collier.hepforge.org/         "COLLIER""
[3]: https://cmake.org/                   "CMake"
