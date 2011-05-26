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
};

#endif
