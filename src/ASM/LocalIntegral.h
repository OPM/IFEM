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

  //! \brief Extracts element solution vectors as \a nsd by \a nen matrices.
  //! \param[in] nsd Number of space dimensions
  //! \param[in] nen Number of element nodes
  //! \param[out] u First element solution vector (deformation)
  //! \param[out] v Second element solution vector (velocity)
  //! \param[out] a Third element solution vector (acceleration)
  //! \param[in] forceCurrent If \e true, return the current values
  //! also if previous values are present
  //!
  //! \details This method is intended for extraction of element-level solution
  //! vectors (deformation, velocity and acceleration for dynamics problems)
  //! on a matrix form, such that equivalent integration point values can be
  //! obtained by multiplying with the basis function value vector, \a fe.N.
  //! If the size of \ref vec is 6 or more and \a forceCurrent is \e false,
  //! the previous values stored in positions {3,4,5} are returned,
  //! otherwise the current values stored in positions {0,1,2} are returned.
  void getSolution(size_t nsd, size_t nen, Matrix* u = nullptr,
                   Matrix* v = nullptr, Matrix* a = nullptr,
                   bool forceCurrent = false) const
  {
    int nvec = vec.size();
    if (nvec > 5 && forceCurrent) nvec = 3;
    int idis = nvec > 5 ? 3 : 0; // index to element displacement vector (u)
    if (nvec == 2 && !forceCurrent) idis = 1; // probably a quasi-static problem
    int ivel = idis + 1; // index to element velocity vector (v)
    int iacc = ivel + 1; // index to element acceleration vector (a)
    if (u && idis < nvec) u->fill(vec[idis],nsd,nen);
    if (v && ivel < nvec) v->fill(vec[ivel],nsd,nen);
    if (a && iacc < nvec) a->fill(vec[iacc],nsd,nen);
  }

  Vectors vec; //!< Element-level primary solution vectors
};

#endif
