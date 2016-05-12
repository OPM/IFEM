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
#include <cstring>

#ifdef HAS_PETSC
#include "petscversion.h"
#ifdef HAVE_MPI
#include "petscsys.h"
#endif
#elif defined(HAVE_MPI)
#include <mpi.h>
#endif

int IFEM::argc;
char** IFEM::argv;
SIMoptions IFEM::cmdOptions;
ControlFIFO IFEM::fifo;
utl::LogStream IFEM::cout(std::cout);


int IFEM::Init (int arg_c, char** arg_v, const char* title)
{
  argc = arg_c;
  argv = arg_v;
  LinAlgInit& linalg = LinAlgInit::Init(argc,argv);
  applyCommandLineOptions(cmdOptions);

  cout.setPIDs(0, linalg.myPid);

  if (linalg.myPid != 0 || argc < 2)
    return linalg.myPid;

  if (title)
  {
    int i, nchar = 13 + strlen(title);
    std::cout <<"\n >>> IFEM "<< title <<" <<<\n ";
    for (i = 0; i < nchar; i++) std::cout <<'=';
    std::cout <<"\n\n Executing command:\n";
#ifdef HAVE_MPI
    int nProc = 1;
#ifdef HAS_PETSC
    MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
    MPI_Comm_size(MPI_COMM_WORLD,&nProc);
#endif
    if (nProc > 1)
      std::cout <<" mpirun -np "<< nProc;
#endif
    for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
    std::cout << std::endl;
  }

  std::cout <<"\n ===== IFEM v"<< IFEM_VERSION_MAJOR <<"."
                               << IFEM_VERSION_MINOR <<"."
                               << IFEM_VERSION_PATCH <<" initialized =====";

  std::cout <<"\n       HDF5 support: ";
#if HAS_HDF5
  std::cout <<"enabled";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n  LR spline support: ";
#if HAS_LRSPLINE
  std::cout <<"enabled";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n     OpenMP support: ";
#if USE_OPENMP
  std::cout <<"enabled";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n        MPI support: ";
#if HAVE_MPI
  std::cout <<"enabled";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n      PETSc support: ";
#if HAS_PETSC
  std::cout <<"enabled (v"<< PETSC_VERSION_MAJOR <<"."
                          << PETSC_VERSION_MINOR <<"."
                          << PETSC_VERSION_SUBMINOR <<")";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n    SuperLU support: ";
#if HAS_SUPERLU
  std::cout <<"enabled (serial)";
#elif HAS_SUPERLU_MT
  std::cout <<"enabled (multi-threaded)";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n       ISTL support: ";
#if HAS_ISTL
  std::cout <<"enabled (v"<< ISTL_VERSION <<")";
#else
  std::cout <<"disabled";
#endif

  std::cout <<"\n        VTF support: ";
#if HAS_VTFAPI == 2
  std::cout <<"enabled (v2)";
#elif HAS_VTFAPI == 1
  std::cout <<"enabled (v1)";
#else
  std::cout <<"disabled";
#endif

  if (cmdOptions.enableController && fifo.open())
    std::cout <<"\n        External controller enabled";

  std::cout << std::endl;

  return 0;
}


void IFEM::applyCommandLineOptions (SIMoptions& opt)
{
  for (int i = 1; i < argc; i++)
    opt.parseOldOptions(argc,argv,i);
}
