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
#include "ControlFIFO.h"
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
#ifdef HAS_UMFPACK
#include <umfpack.h>
#endif


int    IFEM::argc = 0;
char** IFEM::argv = nullptr;
SIMoptions IFEM::cmdOptions;
ControlFIFO* IFEM::fifo = nullptr;
utl::LogStream IFEM::cout(std::cout);
std::shared_ptr<std::ostringstream> IFEM::memoryLog;


int IFEM::Init (int arg_c, char** arg_v, const char* title)
{
  argc = arg_c;
  argv = arg_v;
  LinAlgInit& linalg = LinAlgInit::Init(argc,argv);
  LinAlgInit::increfs();

  bool enableController = false;
  for (int i = 1; i < argc; i++)
    if (!strcasecmp(argv[i],"-controller"))
      enableController = true;
    else
      cmdOptions.parseOldOptions(argc,argv,i);

  cout.setPIDs(0,linalg.myPid);
  memoryLog = std::make_shared<std::ostringstream>();
  cout.addExtraLog(memoryLog);

  if (linalg.myPid != 0 || argc < 2)
    return linalg.myPid;

  if (title)
  {
    int i, nchar = 13 + strlen(title);
    IFEM::cout <<"\n >>> IFEM "<< title <<" <<<\n ";
    for (i = 0; i < nchar; i++) std::cout <<'=';
    IFEM::cout <<"\n\n Executing command:\n";
#ifdef HAVE_MPI
    int nProc = 1;
#ifdef HAS_PETSC
    MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
    MPI_Comm_size(MPI_COMM_WORLD,&nProc);
#endif
    if (nProc > 1)
      IFEM::cout <<" mpirun -np "<< nProc;
#endif
    for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
    IFEM::cout << std::endl;
  }

  IFEM::cout <<"\n ===== IFEM v"<< IFEM_VERSION_MAJOR <<"."
                                << IFEM_VERSION_MINOR <<"."
                                << IFEM_VERSION_PATCH <<" initialized =====";

  IFEM::cout <<"\n       HDF5 support: ";
#if HAS_HDF5
  IFEM::cout <<"enabled";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n  LR spline support: ";
#if HAS_LRSPLINE
  IFEM::cout <<"enabled";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n     OpenMP support: ";
#if USE_OPENMP
  IFEM::cout <<"enabled";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n        MPI support: ";
#ifdef HAVE_MPI
  IFEM::cout <<"enabled";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n      PETSc support: ";
#if HAS_PETSC
  IFEM::cout <<"enabled (v"<< PETSC_VERSION_MAJOR <<"."
                           << PETSC_VERSION_MINOR <<"."
                           << PETSC_VERSION_SUBMINOR <<")";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n    SuperLU support: ";
#if HAS_SUPERLU
  IFEM::cout <<"enabled (serial)";
#elif HAS_SUPERLU_MT
  IFEM::cout <<"enabled (multi-threaded)";
#else
  IFEM::cout <<"disabled";
#endif
  IFEM::cout <<"\n    UMFPack support: ";
#if HAS_UMFPACK
  IFEM::cout <<"enabled (v" << UMFPACK_MAIN_VERSION <<"."
                            << UMFPACK_SUB_VERSION <<"."
                            << UMFPACK_SUBSUB_VERSION << ")";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n       ISTL support: ";
#if HAS_ISTL
  IFEM::cout <<"enabled (v"<< ISTL_VERSION <<")";
#else
  IFEM::cout <<"disabled";
#endif

  IFEM::cout <<"\n        VTF support: ";
#if HAS_VTFAPI == 2
  IFEM::cout <<"enabled (v2)";
#elif HAS_VTFAPI == 1
  IFEM::cout <<"enabled (v1)";
#else
  IFEM::cout <<"disabled";
#endif

  if (enableController)
  {
    fifo = new ControlFIFO();
    if (fifo->open())
      IFEM::cout <<"\n        External controller enabled";
    else
    {
      delete fifo;
      fifo = nullptr;
    }
  }

  IFEM::cout << std::endl;

  return 0;
}


void IFEM::Close ()
{
  delete fifo;
  fifo = nullptr;

  LinAlgInit::decrefs();
}


void IFEM::registerCallback (ControlCallback& cb)
{
  if (fifo) fifo->registerCallback(cb);
}


bool IFEM::pollControllerFifo ()
{
  if (!fifo) return false;

  fifo->poll();
  return true;
}


void IFEM::applyCommandLineOptions (SIMoptions& opt)
{
  for (int i = 1; i < argc; i++)
    opt.parseOldOptions(argc,argv,i);
}
