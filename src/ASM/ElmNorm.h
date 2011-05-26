// $Id$
//==============================================================================
//!
//! \file ElmNorm.h
//!
//! \date Dec 09 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of integrated norm quantities over an element.
//!
//==============================================================================

#ifndef _ELM_NORM_H
#define _ELM_NORM_H

#include "LocalIntegral.h"


/*!
  \brief Class representing integrated norm quantities over an element.
  \details The class is essentially just a double array, but is derived from
  LocalIntegral such that it may be passed as argument to the
  Integrand::evalInt and Integrand::evalBou methods.
*/

class ElmNorm : public LocalIntegral
{
public:
  //! \brief The constructor assigns the internal pointer.
  ElmNorm(double* p) : ptr(p) {}
  //! \brief Empty destructor.
  virtual ~ElmNorm() {}

  //! \brief Indexing operator for assignment.
  double& operator[](size_t i) { return ptr[i]; }
  //! \brief Indexing operator for referencing.
  const double& operator[](size_t i) const { return ptr[i]; }

private:
  double* ptr; //!< Pointer to the actual norm values
};

#endif
