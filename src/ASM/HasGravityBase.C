// $Id$
//==============================================================================
//!
//! \file HasGravityBase.C
//!
//! \date May 27 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class representing FEM integrands with gravity terms.
//!
//==============================================================================

#include "HasGravityBase.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


bool HasGravityBase::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"gravity"))
    return false;
  else if (nsd < 1)
    return true;

  IFEM::cout <<"  Parsing <gravity>\n";
  if (elem->FirstChild())
  {
    char* grav = strdup(elem->FirstChild()->Value());
    char* gval = strtok(grav," ");
    for (unsigned short int d = 0; d < nsd && gval; gval = strtok(nullptr," "))
      gravity[d++] = atof(gval);
    free(grav);
  }

  utl::getAttribute(elem,"x",gravity.x);
  IFEM::cout <<"\tGravitation vector: "<< gravity.x;
  if (nsd >= 2)
  {
    utl::getAttribute(elem,"y",gravity.y);
    IFEM::cout <<" "<< gravity.y;
  }
  if (nsd >= 3)
  {
    utl::getAttribute(elem,"z",gravity.z);
    IFEM::cout <<" "<< gravity.z;
  }
  IFEM::cout << std::endl;

  return true;
}
