// $Id$
//==============================================================================
//!
//! \file AppCommon.h
//!
//! \date Nov 06 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Common helper templates for applications.
//!
//==============================================================================

#ifndef _APP_COMMON_H_
#define _APP_COMMON_H_

#include "XMLInputBase.h"


namespace SIM
{
  //! \brief Base class for input file pre-parsing in applications.
  class AppXMLInputBase : public XMLInputBase
  {
  public:
    //! \brief Default constructor.
    AppXMLInputBase() : dim(3) {}

  protected:
    //! \brief Parses a data section from an XML element.
    virtual bool parse(const TiXmlElement* elem);

  public:
    int dim; //!< Dimensionality of simulation
  };
}

#endif
