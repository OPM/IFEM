// $Id$
//==============================================================================
//!
//! \file LocalIntegral.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing integrated quantities.
//!
//==============================================================================

#ifndef _LOCAL_INTEGRAL_H
#define _LOCAL_INTEGRAL_H

#include "MatVec.h"


/*!
  \brief Abstract base class representing an element level integrated quantity.
*/

class LocalIntegral
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  LocalIntegral() {}

public:
  //! \brief Empty destructor.
  virtual ~LocalIntegral() {}
  //! \brief Virtual destruction method to clean up after numerical integration.
  virtual void destruct() { delete this; }
  //! \brief Returns the LocalIntegral object to assemble into the global one.
  virtual const LocalIntegral* ref() const { return this; }

  Vectors vec; //!< Element-level primary solution vectors
};

#endif
