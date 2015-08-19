//==============================================================================
//!
//! \file ResidualOperators.h
//!
//! \date Aug 19 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete residual operators
//!
//==============================================================================

#ifndef RESIDUALOPERATORS_H_
#define RESIDUALOPERATORS_H_

class FiniteElement;
class Vec3;

#include "MatVec.h"

namespace ResidualOperators
{
  //! \brief Compute a convection term in a residual vector.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] U  Advected field.
  //! \param[in] dUdX Advected field gradient.
  //! \param[in] UC Advecting field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] conservative True to use the conservative formulation.
  void Convection(Vector& EV, const FiniteElement& fe,
                  const Vec3& U, const Matrix& dUdX, const Vec3& UC,
                  double scale, size_t cmp, size_t nf,
                  bool conservative);

  //! \brief Compute a divergence term in a residual vector.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] dUdX Gradient of field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] nf Number of fields in basis.
  //
  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Matrix& dUdX, double scale=1.0,
                  size_t cmp=1, size_t nf=1);

  //! \brief Compute a laplacian term in a residual vector.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] dUdX Current solution gradient.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] stress Whether to add extra stress formulation terms.
  //! \param[in] scmp Starting component.
  void Laplacian(Vector& EV, const FiniteElement& fe,
                 const Matrix& dUdX, double scale=1.0,
                 size_t cmp=1, size_t nf=1,
                 bool stress=false, size_t scmp=0);
}

#endif
