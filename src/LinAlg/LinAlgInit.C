// $Id$
//==============================================================================
//!
//! \file LinAlgInit.C
//!
//! \date May 06 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Initializer for linear algebra packages.
//!
//==============================================================================

#include "LinAlgInit.h"
#include "PETScSupport.h"
#include <iostream>
#if defined(HAVE_MPI) && !defined(HAS_PETSC)
#include <mpi.h>
#endif

LinAlgInit* LinAlgInit::instance = 0;
int         LinAlgInit::refs = 0;


LinAlgInit& LinAlgInit::Init (int argc, char** argv)
{
  if (!instance)
    instance = new LinAlgInit(argc,argv);

  return *instance;
}


LinAlgInit::LinAlgInit (int argc, char** argv)
{
#if defined(HAS_SLEPC)
  SlepcInitialize(&argc,&argv,(char*)0,PETSC_NULL);
#elif defined(HAS_PETSC)
  PetscInitialize(&argc,&argv,(char*)0,PETSC_NULL);
#endif
#ifdef HAVE_MPI
#ifndef HAS_PETSC
  MPI_Init(&argc, &argv);
#endif
  MPI_Comm_rank(MPI_COMM_WORLD,&myPid);
#else
  myPid = 0;
#endif
}


LinAlgInit::~LinAlgInit ()
{
#if defined(HAS_SLEPC)
  SlepcFinalize();
#elif defined(HAS_PETSC)
  PetscFinalize();
#elif defined(HAVE_MPI)
  MPI_Finalize();
#endif
}
