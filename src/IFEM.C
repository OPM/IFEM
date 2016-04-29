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
#ifdef HAS_PETSC
#include "PETScSupport.h"
#endif

int IFEM::argc;
char** IFEM::argv;
SIMoptions IFEM::cmdOptions;
ControlFIFO IFEM::fifo;
utl::LogStream IFEM::cout(std::cout);


int IFEM::Init(int argc_, char** argv_)
{
  argc = argc_;
  argv = argv_;
  LinAlgInit& linalg = LinAlgInit::Init(argc,argv);
  applyCommandLineOptions(cmdOptions);

  cout.setPIDs(0, linalg.myPid);

  if (linalg.myPid != 0 || argc < 2)
    return linalg.myPid;

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

  std::cout <<"\n        MPI support: "<<
#if HAVE_MPI
    "enabled";
#else
    "disabled";
#endif

  std::cout <<"\n      PETSc support: "<<
#if HAS_PETSC
    "enabled (v" << PETSC_VERSION_MAJOR << "."
                 << PETSC_VERSION_MINOR << "."
                 << PETSC_VERSION_SUBMINOR << ")";
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

  if (cmdOptions.enableController && fifo.open())
    std::cout << "\n External controller enabled";

  std::cout << std::endl;

  return 0;
}


void IFEM::applyCommandLineOptions(SIMoptions& opt)
{
  for (int i=1; i < argc; ++i)
    opt.parseOldOptions(argc, argv, i);
}
