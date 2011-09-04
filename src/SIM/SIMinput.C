// $Id$
//==============================================================================
//!
//! \file SIMinput.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for simulators with input parsing functionality.
//!
//==============================================================================

#include "SIMinput.h"
#include "Utilities.h"
#ifdef PARALLEL_PETSC
#include "petscsys.h"
#endif
#include <fstream>


int SIMinput::msgLevel = 2;


SIMinput::SIMinput ()
{
#ifdef PARALLEL_PETSC
  MPI_Comm_rank(PETSC_COMM_WORLD,&myPid);
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
  myPid = 0;
  nProc = 1;
#endif
}


bool SIMinput::read (const char* fileName)
{
  std::ifstream is(fileName);
  if (is)
    std::cout <<"\nReading input file "<< fileName << std::endl;
  else
  {
    std::cerr <<"\n *** SIMinput::read: Failure opening input file "
	      << fileName << std::endl;
    return false;
  }

  char* keyWord = 0;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
		<< keyWord <<"\""<< std::endl;
      return false;
    }

  std::cout <<"\nReading input file succeeded."<< std::endl;
  return true;
}
