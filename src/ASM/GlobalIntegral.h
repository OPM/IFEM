// $Id$
//==============================================================================
//!
//! \file GlobalIntegral.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing integrated quantities.
//!
//==============================================================================

#ifndef _GLOBAL_INTEGRAL_H
#define _GLOBAL_INTEGRAL_H

class LocalIntegral;


/*!
  \brief Abstract base class representing a system level integrated quantity.
*/

class GlobalIntegral
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  GlobalIntegral() {}

public:
  //! \brief Empty destructor.
  virtual ~GlobalIntegral() {}

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  //! \param[in] elmObj The local integral object to add into \a *this.
  //! \param[in] elmId Global number of the element associated with elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId) = 0;
};

#endif
