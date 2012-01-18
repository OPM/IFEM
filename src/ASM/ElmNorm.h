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
#include <cstddef>


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
  //! \param[in] p Pointer to element norm values
  //! \param[in] n Number of norm values
  ElmNorm(double* p, size_t n) : ptr(p), nnv(n) {}
  //! \brief Alternative constructor using the internal buffer \a buf.
  //! \param[in] n Number of norm values
  //!
  //! \details This constructor is used when the element norms are not requested
  //! by the application, but are only used to assembly the global norms.
  //! To avoid the need for a global array of element norms in that case,
  //! an internal array is then used instead.
  ElmNorm(size_t n) : buf(n,0.0), nnv(n) { ptr = &buf.front(); }
  //! \brief Empty destructor.
  virtual ~ElmNorm() {}

  //! \brief Indexing operator for assignment.
  double& operator[](size_t i) { return ptr[i]; }
  //! \brief Indexing operator for referencing.
  const double& operator[](size_t i) const { return ptr[i]; }

  //! \brief Returns the number of norm values.
  size_t size() const { return nnv; }

  //! \brief Returns whether the element norms are stored externally or not.
  bool externalStorage() const { return buf.empty(); }

  //! \brief Virtual destruction method to clean up after numerical integration.
  //! \details Unless the internal buffer \a buf is used, these method only
  //! clears the two solution vector containers (to save memory).
  //! The object itself is then NOT deleted, since it might be used in a second
  //! integration loop over element boundaries, for instance.
  virtual void destruct()
  {
    if (buf.empty())
    {
      vec.clear();
      psol.clear();
    }
    else
      delete this; // The internal buffer has been used, delete the whole thing
  }

private:
  RealArray buf; //!< Internal buffer used when element norms are not requested
  double*   ptr; //!< Pointer to the actual norm values
  size_t    nnv; //!< Number of norm values

public:
  Vectors  psol; //!< Element-level projected solution vectors
};

#endif
