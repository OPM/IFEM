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

#include <cstddef>
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
  ElmNorm(double* p, size_t n) : ptr(p), nnv(n) {}
  //! \brief Empty destructor.
  virtual ~ElmNorm() {}

  //! \brief Indexing operator for assignment.
  double& operator[](size_t i) { return ptr[i]; }
  //! \brief Indexing operator for referencing.
  const double& operator[](size_t i) const { return ptr[i]; }

  //! \brief Returns the number of norm values.
  size_t size() const { return nnv; }

  //! \brief Virtual destruction method to clean up after numerical integration.
  //! \details For this class we only clear the two solution vector containers
  //! (to save memory). We do NOT delete the object itself, since it might be
  //! needed in a second integration loop over element boundaries, for instance.
  virtual void destruct() { vec.clear(); psol.clear(); }

private:
  double* ptr; //!< Pointer to the actual norm values
  size_t  nnv; //!< Number of norm values
};

#endif
