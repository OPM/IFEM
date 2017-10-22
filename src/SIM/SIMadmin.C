// $Id$
//==============================================================================
//!
//! \file SIMadmin.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration base class for FEM simulators.
//!
//==============================================================================

#include "SIMadmin.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <fstream>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif

int SIMadmin::msgLevel = 2;


SIMadmin::SIMadmin (const char* heading) : opt(myOpts)
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


SIMadmin::SIMadmin (SIMadmin& anotherSIM) : opt(anotherSIM.myOpts)
{
  adm = anotherSIM.adm;
  myPid = anotherSIM.myPid;
  nProc = anotherSIM.nProc;
}


void SIMadmin::printHeading (int& subStep) const
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


bool SIMadmin::read (const char* fileName)
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


bool SIMadmin::readFlat (const char* fileName)
{
  std::ifstream is(fileName);
  if (!is)
  {
    std::cerr <<"\n *** SIMadmin::read: Failure opening input file \""
              << fileName <<"\"."<< std::endl;
    return false;
  }

  IFEM::cout <<"\nReading input file "<< fileName << std::endl;

  char* keyWord = nullptr;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                << keyWord <<"\"."<< std::endl;
      return false;
    }

  IFEM::cout <<"\nReading input file succeeded."<< std::endl;

  return true;
}


bool SIMadmin::parse (char*, std::istream&)
{
  std::cerr <<" *** SIMadmin::parse(char*,std::istream&):"
            <<" The flat file format is depreciated.\n"
            <<"     Use the XML format instead."<< std::endl;
  return false;
}


bool SIMadmin::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"messageLevel");
  if (value) msgLevel = atoi(value);
  return true;
}
