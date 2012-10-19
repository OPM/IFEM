//==============================================================================
//!
//! \file IFEM.C
//!
//! \date Aug 08 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Initialization of the IFEM library.
//!
//==============================================================================

#include "IFEM.h"
#include "LinAlgInit.h"

#include <iostream>

SIMoptions IFEM_cmdOptions;
int IFEM_argc;
char** IFEM_argv;

bool InitIFEM(int argc, char** argv, int myId)
{
  IFEM_argc = argc;
  IFEM_argv = argv;
  LinAlgInit::Init(argc,argv);
  for (int i=1; i < argc; ++i)
    IFEM_cmdOptions.parseOldOptions(argc, argv, i);

  if (myId != 0)
    return true;

  std::cout << "\n===== IFEM v" << IFEM_VERSION_MAJOR << "." 
                                << IFEM_VERSION_MINOR << "." 
                                << IFEM_VERSION_PATCH << " initialized =====" << std::endl;

  std::cout << "       HDF5 support: " <<
#if HAS_HDF5
    "enabled\n";
#else
    "disabled\n";
#endif

  std::cout << "  LR spline support: " <<
#if HAS_LRSPLINE
    "enabled\n";
#else
    "disabled\n";
#endif

  std::cout << "     OpenMP support: " <<
#if USE_OPENMP
    "enabled\n";
#else
    "disabled\n";
#endif

  std::cout << "      PETSc support: " <<
#if PARALLEL_PETSC
    "enabled (MPI)\n";
#elif HAS_PETSC
    "enabled\n";
#else
    "disabled\n";
#endif

  std::cout << "    SuperLU support: " << 
#if HAS_SUPERLU
    "enabled (serial)\n";
#elif HAS_SUPERLU_MT
    "enabled (multi-threaded)\n";
#else
    "disabled\n";
#endif

  std::cout << "        VTF support: " << 
#if HAS_VTFAPI == 2
    "enabled (v2)\n";
#elif HAS_VTFAPI == 1
    "enabled (v1)\n";
#else
    "disabled\n";
#endif

  return true;
}
