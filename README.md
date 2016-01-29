# IFEM project


## Introduction

IFEM is an object-oriented toolbox for implementing isogeometric finite element
solvers for linear and nonlinear partial differential equations.
The toolbox is developed through the ICADA project at SINTEF ICT, Trondheim.
The purpose of this code is to serve as a common base for isogeometric
PDE-simulators, using splines and NURBS as basis functions in the finite element
formulations. The toolbox contains methods for doing linear and non-linear,
stationary and dynamic time-domain analyses, as well as eigenvalue analyses.

## Module overview

The simulation toolbox is organized into a set of modules,
organized as class hierarchies implemented in the C++ language.
The top-level driver is organized in the class SIMbase and its sub-classes.
These classes have methods for reading model data from an input file,
to assemble the linearized algebraic system of equations resulting from the
finite element discretization, and to solve that system.
It also contains methods for writing a VTF-file with results.
Problem-specific drivers and main programs are found in the sub-folder Apps.

The core of the finite element implementation is contained in the class
ASMbase and its sub-classes, which have methods for evaluating the element
matrices involved and assembling them into the system matrices.
There is typically one ASMbase object for each spline patch in the model.

The physical problem-dependent data and methods is accessed via an abstract
interface class, Integrand, through which the application programmer can
implement the weak form of the underlying finite element problem.

The actual splines evaluation is performed through the GoTools library, which is
not part of the current package. In addition, this code depends on the ARPACK,
LAPACK and BLAS libraries. Optionally, the SAMG algebraic multi-grid solver and
the SuperLU direct solver (public domain http://crd.lbl.gov/~xiaoye/SuperLU)
may be included. The Ceetron VTFAPI library may also be included for direct
export of simulations results to GLview VTF-files.

### Getting all dependencies

A number of things need to be set up properly to compile the IFEM library.
First, we will add inhouse dependencies by a secondary repository.

1. Add the IFEM repository at https://launchpad.net/~ifem/ (follow the instructions on site)
2. Install development tools and compilers by typing

    `sudo apt-get install cmake g++ gfortran`

3. Install official libraries

    `sudo apt-get install python-dev libnewmat10ldbl libboost-dev libblas-dev liblapack-dev libarpack2-dev libsuperlu3-dev`

4. Install inhouse libraries

    `sudo apt-get install libgotools-compositemodel-dev libgotools-core-dev libgotools-igeslib-dev libgotools-implicitization-dev libgotools-intersections-dev libgotools-isogeometricmodel-dev libgotools-qualitymodule-dev libgotools-qualitymodule1 libgotools-parametrization-dev libgotools-topology-dev libgotools-trivariate-dev libgotools-trivariatemodel-dev libttl-dev libsisl-dev`

5. **[optional]** Install petsc by the official webpage (http://www.mcs.anl.gov/petsc/download/)
6. **[optional]** Install LR-splines by typing

    `sudo apt-get install liblrspline1-dev`

7. **[optional]** Install VTF writer. This is proprietary software and cannot be shared openly. E-mail Trond.Kvamsdal@sintef.no and ask for them.


### Getting the code

This is done by first navigating to the folder in which you want IFEM installed and typing

    git clone https://github.com/OPM/IFEM


### Compiling the code

To compile, first navigate to the root catalogue of IFEM, here denoted by `<IFEM root>`.

1. `cd <IFEM root>`
2. `mkdir Debug`
3. `cd Debug`
4. **[optional]** specify which submodules you have available in `<IFEM root>/cmake/Modules/IFEMoptions.cmake`
5. `cmake -DCMAKE_BUILD_TYPE=Debug ..`
6. `make `

this will compile the library.
Change all instances of `Debug` with `Release` to drop debug-symbols, but get faster running code.


### Testing the code

IFEM is using cmake test system. To compile run all regression- and unit-tests, navigate to your build
folder (i.e. `<IFEM root>/Debug`) and type

    make check


