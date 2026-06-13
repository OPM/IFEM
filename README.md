## Introduction

**IFEM** is an object-oriented toolbox for implementing isogeometric finite
element (FE) solvers for linear and nonlinear partial differential equations.
The main toolbox together with structural mechanics applications was developed
through the ICADA project at SINTEF Digital in collaboration with NTNU.
Important extensions to fluid mechanics and coupled problems was done
during the NOWITECH and FSI-WT projects.

The purpose of this project is to serve as a common base for isogeometric
PDE-simulators, using splines and NURBS as basis functions in the FE
formulations. The toolbox contains methods for doing linear and non-linear,
stationary and dynamic time-domain analyses, as well as eigenvalue analyses.
**IFEM** is parallelized using the [PETSc](https://petsc.org) library.
Adaptive refinement is enabled using [LR B-splines](https://lrbsplines.com)
based on residual or recovery methods for a posteriori error estimation.

## Module overview

The **IFEM** toolbox is organized into a set of modules,
organized as class hierarchies implemented in the C++ language.
The top-level driver is organized in the class `SIMbase` and its sub-classes.
These classes have methods for reading model data from an input file,
to assemble the linearized algebraic system of equations resulting from the
finite element discretization, and to solve that system.
It also contains methods for exporting the simulation results to HDF5 files
for post-processing in external modules.

The core of the FE implementation is contained in the class `ASMbase` and
its sub-classes, which have methods for evaluating the element matrices
involved and assembling them into the system matrices.
There is typically one `ASMbase` object for each spline patch in the model.

The physical problem-dependent data and methods are accessed via an abstract
interface class, `Integrand`, through which the application programmer can
implement the weak form of the underlying finite element problem.

Refer to the [github pages](https://opm.github.io/IFEM)
for detailed documentation on the various class hierarchies in **IFEM**.

The actual splines evaluation in **IFEM** is performed through the
[GoTools](https://github.com/SINTEF-Geometry/GoTools) library,
which is not part of this project.
In addition, the toolbox depends on the ARPACK, LAPACK and BLAS libraries
for doing low-level linear algebra operations.
[SuperLU](https://portal.nersc.gov/project/sparse/superlu/)
is used for direct solution of the linear equation systems by default.
Other equation solvers may be used by including some of the optional packages
listed below.

## Getting all dependencies

A number of things need to be set up properly to build the IFEM library:

1. Add the IFEM personal package archive

       sudo apt-add-repository ppa:ifem/ppa
       sudo apt-get update

   See also https://launchpad.net/~ifem/.

2. Install development tools and compilers by typing

       sudo apt-get install cmake g++ gfortran

3. Install required official libraries

       sudo apt-get install libboost-dev libarpack2-dev libsuperlu-dev catch2 libtinyxml2-dev

4. Install required in-house GoTools libraries by typing

       sudo apt-get install libgotools-core-dev libgotools-trivariate-dev

   There are a number of other GoTools libraries available also, but
   the two above are the only ones required to build IFEM simulators.

5. **[optional]** Install PETSc from the official web page
   https://petsc.org/release/install/download

6. **[optional]** Install UMFPACK direct solver

       sudo apt-get install libsuitesparse-dev

7. **[optional]** Install LR-splines by typing

       sudo apt-get install liblrspline1-dev

8. **[optional]** Install support for HDF5 output typing (one or both commands,
   the second one is needed only if you want to build parallel applications)

       sudo apt-get install libhdf5-serial-dev
       sudo apt-get install libhdf5-openmpi-dev

9. **[optional]** Install support for [de]serialization of simulator state.
   This package is needed if you want to equip your simulator with restart
   capabilities. Install it by typing

       sudo apt-get install libcereal-dev

10. **[optional]** Install the Ceetron VTFAPI library.
   This is proprietary software that cannot be shared openly.  
   Send email to Trond.Kvamsdal@sintef.no and ask for them.
   This library enables direct export of simulations results to VTF-files
   for visualization in the proprietary GLview Inova software.

## Getting the code

Navigate to the folder in which you want the IFEM source installed and type

    git clone https://github.com/OPM/IFEM

## Compiling the code

Navigate to the root folder of IFEM source, here denoted by `<IFEM root>`. Then

    mkdir Debug
    cd Debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug [<IFEM-options>]
    make

where `<IFEM-options>` is an optional list of sub-modules
you may choose to activate or deactivate in your build.
See the file `<IFEM root>/cmake/Modules/IFEMOptions.cmake` for a complete list
of available options and their default setting.

This will build the libraries which can be found in the `Debug/lib` sub-folder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

## Testing the code

IFEM uses the cmake test system (CTest).
To build and run all regression- and unit-tests, navigate to your build
folder (i.e., `<IFEM root>/Debug`) and type

    make check
