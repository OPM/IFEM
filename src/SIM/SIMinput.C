// $Id: SIMinput.C,v 1.3 2010-10-10 11:20:54 kmo Exp $
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
  char* keyWord = 0;
  std::cout <<"\nReading input file "<< fileName << std::endl;
  std::ifstream is(fileName);
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
