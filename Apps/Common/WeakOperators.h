//==============================================================================
//!
//! \file WeakOperators.h
//!
//! \date Jul 22 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete operators.
//!
//==============================================================================

#ifndef WEAKOPERATORS_H_
#define WEAKOPERATORS_H_

class FiniteElement;
class Vec3;

#include "MatVec.h"

namespace WeakOperators
{
  //! \brief Compute an advection term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] AC Advecting field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scmp Starting component.
  void Advection(Matrix& EM, const FiniteElement& fe,
                 const Vec3& AC, double scale=1.0,
                 size_t cmp=1, size_t nf=1, size_t scmp=0);

  //! \brief Compute an (nonlinear) convection term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] U  Advecting field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] conservative True to use the conservative formulation.
  void Convection(Matrix& EM, const FiniteElement& fe,
                  const Vec3& U, const Matrix& dUdX, double scale,
                  size_t cmp, size_t nf, bool conservative);

  //! \brief Compute a divergence term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scale Scaling factor for contribution.
  void Divergence(Matrix& EM, const FiniteElement& fe,
                  size_t nf, double scale=1.0);

  //! \brief Compute a divergence term.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] D Divergence of field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] nf Number of fields in basis.
  //
  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Vec3& D, double scale=1.0,
                  size_t cmp=1, size_t nf=1);

  //! \brief Compute a divergence/gradient term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scale Scaling factor for contribution.
  void PressureDiv(Matrix& EM, const FiniteElement& fe,
                   size_t nf, double scale=1.0);

  //! \brief Compute a gradient term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scale Scaling factor for contribution.
  void Gradient(Matrix& EM, const FiniteElement& fe,
                size_t nf, double scale=1.0);

  //! \brief Compute a gradient term.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] nf Number of fields in basis.
  void Gradient(Vector& EV, const FiniteElement& fe,
                double scale=1.0, size_t nf=1);

  //! \brief Compute a laplacian.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] stress Whether to add extra stress formulation terms.
  //! \param[in] scmp Starting component.
  void Laplacian(Matrix& EM, const FiniteElement& fe,
                 double scale=1.0, size_t cmp=1, size_t nf=1,
                 bool stress=false, size_t scmp=0);

  //! \brief Compute a heteregenous coefficient laplacian.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[out] K The coefficient matrix.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Scaling factor for contribution.
  void LaplacianCoeff(Matrix& EM, const Matrix& K, const FiniteElement& fe,
                      double scale=1.0);

  //! \brief Compute a mass term.
  //! \param[out] EM The element matrix to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scmp Starting component.
  void Mass(Matrix& EM, const FiniteElement& fe,
            double scale=1.0, size_t cmp=1, size_t nf=1, size_t scmp=0);

  //! \brief Compute a source term.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scmp Starting component.
  void Source(Vector& EV, const FiniteElement& fe,
              double scale=1.0, size_t cmp=1, size_t nf=1, size_t scmp=0);

  //! \brief Compute a vector-source term.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] scale Vector with contributions.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] scmp Starting component.
  void Source(Vector& EV, const FiniteElement& fe,
              const Vec3& f, double scale=1.0, size_t cmp=1, size_t nf=1, size_t scmp=0);
}

#endif
