// $Id: LinAlgInit.C,v 1.2 2011-01-02 15:50:35 kmo Exp $
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
#ifdef HAS_PETSC
#include "slepceps.h"
#include "petscksp.h"
#endif


LinAlgInit::LinAlgInit (int argc, char** argv)
{
#ifdef HAS_PETSC
  SlepcInitialize(&argc,&argv,(char*)0,PETSC_NULL);
  //PetscInitialize(&argc,&argv,(char*)0,PETSC_NULL);
#endif
#ifdef PARALLEL_PETSC
  MPI_Comm_rank(PETSC_COMM_WORLD,&myPid);
#else
  myPid = 0;
#endif
}


LinAlgInit::~LinAlgInit ()
{
#ifdef HAS_PETSC
  SlepcFinalize();
  //PetscFinalize();
#endif
}
