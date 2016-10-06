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

class Vec3;

#include "FiniteElement.h"
#include "MatVec.h"

namespace WeakOperators
{
  //! \brief Compute an advection term.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] AC Advecting field
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis to use
  void Advection(Matrix& EM, const FiniteElement& fe,
                 const Vec3& AC, double scale=1.0, int basis=1);

  //! \brief Compute a (nonlinear) convection term.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] U  Advecting field
  //! \param[in] conservative True to use the conservative formulation
  //! \param[in] basis Basis to use
  template<class T>
  void Convection(Matrix& EM, const FiniteElement& fe,
                  const Vec3& U, const T& dUdX, double scale,
                  bool conservative=false, int basis=1)
  {
    size_t cmp = EM.rows() / fe.basis(basis).size();
    if (conservative) {
      Advection(EM, fe, U, -scale, basis);
      for (size_t i = 1;i <= fe.basis(basis).size();i++)
        for (size_t j = 1;j <= fe.basis(basis).size();j++) {
          for (size_t k = 1;k <= cmp;k++) {
            for (size_t l = 1;l <= cmp;l++)
              EM((j-1)*cmp+l,(i-1)*cmp+k) -= scale*U[l-1]*fe.basis(basis)(i)*fe.grad(basis)(j,k)*fe.detJxW;
          }
        }
    }
    else {
      Advection(EM, fe, U, scale, basis);
      for (size_t i = 1;i <= fe.basis(basis).size();i++)
        for (size_t j = 1;j <= fe.basis(basis).size();j++) {
          for (size_t k = 1;k <= cmp;k++) {
            for (size_t l = 1;l <= cmp;l++)
              EM((j-1)*cmp+l,(i-1)*cmp+k) += scale*dUdX(l,k)*fe.basis(basis)(j)*fe.basis(basis)(i)*fe.detJxW;
          }
        }
    }
  }

  //! \brief Compute a divergence term.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis for field
  //! \param[in] tbasis Test function basis
  void Divergence(Matrix& EM, const FiniteElement& fe,
                  double scale=1.0, int basis=1, int tbasis=1);

  //! \brief Compute a divergence term.
  //! \param[out] EV The element vector to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] D Divergence of field
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Test function basis
  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Vec3& D, double scale=1.0, int basis=1);

  //! \brief Compute a gradient term for a (potentially mixed) vector/scalar field.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis for field
  //! \param[in] tbasis Test function basis
  void Gradient(Matrix& EM, const FiniteElement& fe,
                double scale=1.0, int basis=1, int tbasis=1);

  //! \brief Compute a gradient term.
  //! \param[out] EV The element vector to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] tbasis Test function basis
  void Gradient(Vector& EV, const FiniteElement& fe,
                double scale=1.0, int basis=1);

  //! \brief Compute a laplacian.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] stress Whether to add extra stress formulation terms
  //! \param[in] basis Basis to use
  void Laplacian(Matrix& EM, const FiniteElement& fe,
                 double scale=1.0, bool stress=false, int basis=1);

  //! \brief Compute a heteregenous coefficient laplacian.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[out] K The coefficient matrix
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  void LaplacianCoeff(Matrix& EM, const Matrix& K, const FiniteElement& fe,
                      double scale=1.0, int basis=1);

  //! \brief Compute a mass term.
  //! \param[out] EM The element matrix to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis to use
  void Mass(Matrix& EM, const FiniteElement& fe,
            double scale=1.0, int basis=1);

  //! \brief Compute a source term.
  //! \param[out] EV The element vector to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis to use
  //! \param[in] cmp Component to add (0 for all)
  void Source(Vector& EV, const FiniteElement& fe,
              double scale=1.0, int cmp=1, int basis=1);

  //! \brief Compute a vector-source term.
  //! \param[out] EV The element vector to add contribution to
  //! \param[in] fe The finite element to evaluate for
  //! \param[in] scale Vector with contributions
  //! \param[in] scale Scaling factor for contribution
  //! \param[in] basis Basis to use
  void Source(Vector& EV, const FiniteElement& fe,
              const Vec3& f, double scale=1.0, int basis=1);
}

#endif
