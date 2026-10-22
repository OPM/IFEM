// $Id$
//==============================================================================
//!
//! \file SIMargsBase.C
//!
//! \date Jul 15 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for pre-parsing of XML input files for simulators.
//!
//==============================================================================

#include "SIMargsBase.h"
#include "SIMoptions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <cstring>
#include <cstdlib>
#include <cctype>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif


bool SIMargsBase::parseArg (const char* argv)
{
  if (!strncmp(argv,"-adap",5))
    adap = strlen(argv) > 5 ? 1+atoi(argv+5) : 1;
  else if (!strncmp(argv,"-noadap",7))
    adap = 0;
  else if (!strcmp(argv,"-1D"))
    dim = 1;
  else if (!strcmp(argv,"-2D"))
    dim = 2;
  else if (!strcmp(argv,"-3D"))
    dim = 3;
  else
    return false;

  return true;
}


bool SIMargsBase::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),context))
    if (utl::getAttribute(elem,"adaptive",adap,false))
      if (adap == '0' || tolower(adap) == 'n')
        adap = 0;

  if (strcasecmp(elem->Value(),"geometry"))
    return true;

  if (!utl::getAttribute(elem,"dimension",dim))
    utl::getAttribute(elem,"dim",dim);

  return IFEM::getOptions().parseDiscretizationTag(elem);
}


/*!
  \details A null pointer is returned if no command-line arguments are specified
  or they all match a defined command-line option. The first argument which does
  not match a defined command-line option is assumed to be the input file name.
  \note If a parsing failure is detected, this method will abort the program
  by invoking exit(1).
*/

char* SIMargsBase::preParse (int argc, char** argv,
                             std::function<bool(int&)> external_parser)
{
  SIMoptions ignored;
  char* input_file = nullptr;
  for (int i = 1; i < argc; i++)
    if (argv[i] == input_file)
      ; // ignore the input file on the second pass
    else if (external_parser(i) || this->parseArg(argv[i]))
      ; // this option was parsed
    else if (ignored.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (input_file)
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;
    else if (strcasestr(input_file = argv[i],".xinp"))
    {
      if (this->readXML(input_file,false))
        i = 0; // start over and let command-line options override input file
      else
        exit(1); // pre-parse failure, terminate immediately
    }

  return input_file;
}
