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
#include "IFEM.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "SystemMatrix.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <fstream>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif

int SIMinput::msgLevel = 2;


SIMinput::SIMinput (const char* heading) : opt(myOpts)
#ifdef HAS_PETSC
  , adm(PETSC_COMM_WORLD)
#elif defined(HAVE_MPI)
  , adm(true)
#endif
{
  myPid = adm.getProcId();
  nProc = adm.getNoProcs();

  myOpts = IFEM::getOptions(); // Initialize options from command-line arguments

  if (heading) myHeading = heading;
}


SIMinput::SIMinput (SIMinput& anotherSIM) : opt(anotherSIM.myOpts)
{
  adm = anotherSIM.adm;
  myPid = anotherSIM.myPid;
  nProc = anotherSIM.nProc;
}


void SIMinput::printHeading (int& subStep) const
{
  if (myHeading.empty())
    return;

  size_t n = myHeading.find_last_of('\n');
  if (n+1 < myHeading.size()) n = myHeading.size()-n;
  IFEM::cout <<"\n"<< ++subStep <<". "<< myHeading <<"\n";
  for (size_t i = 0; i < 3+n && i < 3+myHeading.size(); i++) IFEM::cout <<'=';
  if (subStep > 9) IFEM::cout <<'=';
  IFEM::cout << std::endl;
}


bool SIMinput::read (const char* fileName)
{
  static int substep = 0;
  this->printHeading(substep);

  bool result;
  if (strcasestr(fileName,".xinp"))
    result = this->readXML(fileName);
  else
    result = this->readFlat(fileName);

  // Let command-line options override settings on the input file
  IFEM::applyCommandLineOptions(opt);

  return result;
}


bool SIMinput::readFlat (const char* fileName)
{
  std::ifstream is(fileName);
  if (!is)
  {
    std::cerr <<"\n *** SIMinput::read: Failure opening input file "
              << fileName << std::endl;
    return false;
  }

  IFEM::cout <<"\nReading input file "<< fileName << std::endl;

  char* keyWord = 0;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                << keyWord <<"\""<< std::endl;
      return false;
    }

  IFEM::cout <<"\nReading input file succeeded."<< std::endl;

  return true;
}


bool SIMinput::parse (char*, std::istream&)
{
  std::cerr <<" *** SIMinput::parse(char*,std::istream&):"
            <<" The flat file format is depreciated.\n"
            <<"     Use the XML format instead."<< std::endl;
  return false;
}
