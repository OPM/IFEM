// $Id$
//==============================================================================
//!
//! \file AppCommon.C
//!
//! \date Nov 06 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Common helper templates for applications.
//!
//==============================================================================

#include "AppCommon.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


bool SIM::AppXMLInputBase::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"geometry"))
    return true;

  if (!utl::getAttribute(elem,"dimension",dim))
    utl::getAttribute(elem,"dim",dim);

  std::string type;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"patchfile"))
      if (utl::getAttribute(child,"type",type) && type == "lrspline")
        IFEM::getOptions().discretization = ASM::LRSpline;

  return true;
}
