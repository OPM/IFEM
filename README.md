# IFEM project

## Introduction

IFEM is an object-oriented toolbox for implementing isogeometric finite element
solvers for linear and nonlinear partial differential equations.
The main toolbox together with structural mechanics applications was developed
through the ICADA project at SINTEF Digital in collaboration with NTNU.
The important extension to fluid mechanics and coupled problems was done
during the NOWITECH and FSI-WT projects.

The purpose of this code is to serve as a common base for isogeometric
PDE-simulators, using splines and NURBS as basis functions in the finite element
formulations. The toolbox contains methods for doing linear and non-linear,
stationary and dynamic time-domain analyses, as well as eigenvalue analyses.
IFEM is parallelized using the PETSc library and enables adaptive refinement
using LR B-splines based on residual or recovery methods for a posteriori error estimation.

## Module overview

The simulation toolbox is organized into a set of modules,
organized as class hierarchies implemented in the C++ language.
The top-level driver is organized in the class `SIMbase` and its sub-classes.
These classes have methods for reading model data from an input file,
to assemble the linearized algebraic system of equations resulting from the
finite element discretization, and to solve that system.
It also contains methods for writing a VTF-file with results.

The core of the finite element implementation is contained in the class
`ASMbase` and its sub-classes, which have methods for evaluating the element
matrices involved and assembling them into the system matrices.
There is typically one `ASMbase` object for each spline patch in the model.

The physical problem-dependent data and methods are accessed via an abstract
interface class, `Integrand`, through which the application programmer can
implement the weak form of the underlying finite element problem.

The actual splines evaluation is performed through the GoTools library, which is
not part of the current package. In addition, this code depends on the ARPACK,
LAPACK and BLAS libraries. SuperLU (http://crd.lbl.gov/~xiaoye/SuperLU)
is used for direct solution of the linear equation systems.
The Ceetron VTFAPI library may also be included for direct
export of simulations results to GLview VTF-files.
Several other optional packages may also be included, as indicated below.

### Getting all dependencies

A number of things need to be set up properly to compile the IFEM library.
First, we will add in-house dependencies by a secondary repository.

1. Add the IFEM repository at https://launchpad.net/~ifem/
   (follow the instructions on that site)

2. Install development tools and compilers by typing

    `sudo apt-get install cmake g++ gfortran`

3. Install required official libraries

    `sudo apt-get install libboost-dev libarpack2-dev libsuperlu-dev libgtest-dev`

4. Install required in-house GoTools libraries by typing

    `sudo apt-get install libgotools-core-dev libgotools-trivariate-dev`

   There are a number of other GoTools libraries available also, but
   the two above are the only ones required to build IFEM simulators.

5. **[optional]** Install PETSc from the official web page
   http://www.mcs.anl.gov/petsc/download/

6. **[optional]** Install UMFPACK direct solver

    `sudo apt-get install libsuitesparse-dev`

7. **[optional]** Install LR-splines by typing

    `sudo apt-get install liblrspline1-dev`

8. **[optional]** Install support for HDF5 output typing (one or both commands,
   the second one is needed only if you want to build parallel applications)

    `sudo apt-get install libhdf5-serial-dev`  
    `sudo apt-get install libhdf5-openmpi-dev`

9. **[optional]** Install support for [de]serialization of simulator state.
   This package is needed if you want to equip your simulator with restart
   capabilities. Install it by typing

    `sudo apt-get install libcereal-dev`

    Note: it is only available for Ubuntu 16.04 and later.

10. **[optional]** Install the VTF writer.
   This is proprietary software that cannot be shared openly.  
   Send email to Trond.Kvamsdal@sintef.no and ask for them.

### Getting the code

Navigate to the folder in which you want the IFEM source installed and type

    git clone https://github.com/OPM/IFEM

### Compiling the code

Navigate to the root folder of IFEM source, here denoted by `<IFEM root>`. Then

1. `mkdir Debug`
2. `cd Debug`
3. `cmake .. -DCMAKE_BUILD_TYPE=Debug [<IFEM-options>]`
4. `make `

where `<IFEM-options>` is an optional list of sub-modules
you may choose to activate or deactivate in your build.
See the file `<IFEM root>/cmake/Modules/IFEMOptions.cmake` for a complete list
of available options and their default setting.

This will compile the libraries which can be found in the `Debug/lib`sub-folder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

### Testing the code

IFEM uses the cmake test system.
To compile and run all regression- and unit-tests, navigate to your build
folder (i.e. `<IFEM root>/Debug`) and type

    make check
