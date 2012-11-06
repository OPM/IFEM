// $Id$
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
int        IFEM_argc;
char**     IFEM_argv;


int InitIFEM (int argc, char** argv)
{
  IFEM_argc = argc;
  IFEM_argv = argv;
  int  myId = LinAlgInit::Init(argc,argv).myPid;

  for (int i=1; i < argc; ++i)
    IFEM_cmdOptions.parseOldOptions(argc, argv, i);

  if (myId != 0)
    return myId;

  std::cout <<"\n ===== IFEM v"<< IFEM_VERSION_MAJOR <<"."
                               << IFEM_VERSION_MINOR <<"."
                               << IFEM_VERSION_PATCH <<" initialized =====";

  std::cout <<"\n       HDF5 support: "<<
#if HAS_HDF5
    "enabled";
#else
    "disabled";
#endif

  std::cout <<"\n  LR spline support: "<<
#if HAS_LRSPLINE
    "enabled";
#else
    "disabled";
#endif

  std::cout <<"\n     OpenMP support: "<<
#if USE_OPENMP
    "enabled";
#else
    "disabled";
#endif

  std::cout <<"\n      PETSc support: "<<
#if PARALLEL_PETSC
    "enabled (MPI)";
#elif HAS_PETSC
    "enabled";
#else
    "disabled";
#endif

  std::cout <<"\n    SuperLU support: "<<
#if HAS_SUPERLU
    "enabled (serial)";
#elif HAS_SUPERLU_MT
    "enabled (multi-threaded)";
#else
    "disabled";
#endif

  std::cout <<"\n        VTF support: "<<
#if HAS_VTFAPI == 2
    "enabled (v2)";
#elif HAS_VTFAPI == 1
    "enabled (v1)";
#else
    "disabled";
#endif

  std::cout << std::endl;
  return 0;
}

//TODO (kmo): Remove this and update all main programs accordingly...
bool InitIFEM(int argc, char** argv, int) { return InitIFEM(argc,argv); }
