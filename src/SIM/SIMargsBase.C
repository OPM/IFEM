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
#include <cctype>


bool SIMargsBase::parseArg (const char* argv)
{
  if (!strncmp(argv,"-adap",5) && !isdigit(argv[5]))
    adap = true;
  else if (!strncmp(argv,"-noadap",7))
    adap = false;
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
    utl::getAttribute(elem,"adaptive",adap);

  if (strcasecmp(elem->Value(),"geometry"))
    return true;

  if (!utl::getAttribute(elem,"dimension",dim))
    utl::getAttribute(elem,"dim",dim);

  return IFEM::getOptions().parseDiscretizationTag(elem);
}
