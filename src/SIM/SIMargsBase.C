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
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <cstring>
#include <cstdlib>
#include <cctype>


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


bool SIMargsBase::parse (const TiXmlElement* elem)
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
