// $Id$
//==============================================================================
//!
//! \file IFEM_math.h
//!
//! \date Oct 20 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various math utility methods.
//!
//==============================================================================

#ifndef UTL_MATH_H
#define UTL_MATH_H

#include <cmath>
#include <cstddef>
#include <vector>


namespace utl
{
  //! \brief Returns the number of monomials in Pascal's triangle.
  //! \param[in] p Polynomial order (>= 0)
  //! \param[in] nsd Number of spatial dimensions (2 or 3)
  size_t Pascal(int p, unsigned short int nsd);
  //! \brief Evaluates the monomials of Pascal's triangle in 2D for order \a p.
  void Pascal(int p, Real x, Real y, std::vector<Real>& phi);
  //! \brief Evaluates the monomials of Pascal's triangle in 3D for order \a p.
  void Pascal(int p, Real x, Real y, Real z, std::vector<Real>& phi);

  //! \brief Evaluates the positive ramp function \f$(x+|x|)/2\f$.
  inline Real Pos(Real x) { return x > Real(0) ? x : Real(0); }
  //! \brief Evaluates the negative ramp function \f$(x-|x|)/2\f$.
  inline Real Neg(Real x) { return x < Real(0) ? x : Real(0); }
  //! \brief Converts from degrees to radians.
  inline Real Rad(Real x) { return x*M_PI/Real(180); }
}

#endif
