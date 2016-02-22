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

class Vec3;

#include "FiniteElement.h"
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
  template<class T>
  void Convection(Vector& EV, const FiniteElement& fe,
                  const Vec3& U, const T& dUdX, const Vec3& UC,
                  double scale, size_t cmp, size_t nf,
                  bool conservative)
  {
    if (conservative) {
      for (size_t i = 1;i <= fe.N.size();i++)
        for (size_t k = 1;k <= cmp;k++)
          for (size_t l = 1;l <= cmp;l++)
            // Convection
            EV((i-1)*nf + k) += scale*U[k-1]*UC[l-1]*fe.dNdX(i,l)*fe.detJxW;
    }
    else {
      for (size_t k = 1;k <= cmp;k++) {
        // Convection
        double conv = 0.0;
        for (size_t l = 1;l <= cmp;l++)
          conv += UC[l-1]*dUdX(k,l);
        conv *= scale*fe.detJxW;

        for (size_t i = 1;i <= fe.N.size();i++)
          EV((i-1)*nf + k) -= conv*fe.N(i);
      }
    }
  }

  //! \brief Compute a divergence term in a residual vector.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] dUdX Gradient of field.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] nf Number of fields in basis.
  template<class T>
  void Divergence(Vector& EV, const FiniteElement& fe,
                  const T& dUdX, double scale=1.0,
                  size_t cmp=1, size_t nf=1, size_t basis=1, size_t ofs=0)
  {
    for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
      double div=0.0;
      for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
        div += dUdX(k,k);
      EV(ofs+(i-1)*nf+cmp) += scale*div*fe.basis(basis)(i)*fe.detJxW;
    }
  }

  //! \brief Compute a laplacian term in a residual vector.
  //! \param[out] EV The element vector to add contribution to.
  //! \param[in] fe The finite element to evaluate for.
  //! \param[in] dUdX Current solution gradient.
  //! \param[in] scale Scaling factor for contribution.
  //! \param[in] cmp Number of components to add.
  //! \param[in] nf Number of fields in basis.
  //! \param[in] stress Whether to add extra stress formulation terms.
  //! \param[in] scmp Starting component.
  template<class T>
  void Laplacian(Vector& EV, const FiniteElement& fe,
                 const T& dUdX, double scale=1.0,
                 size_t cmp=1, size_t nf=1,
                 bool stress=false, size_t scmp=0)
  {
    for (size_t i = 1;i <= fe.N.size();i++) {
      for (size_t k = 1;k <= cmp;k++) {
        double diff = 0.0;
        for (size_t l = 1;l <= cmp;l++)
          diff += dUdX(k,l)*fe.dNdX(i,l);
        diff *= scale*fe.detJxW;

        // Add negative residual to rhs of momentum equation
        EV((i-1)*nf + k+scmp) -= diff;
      }
    }

    // Use stress formulation
    if (stress) {
      for (size_t i = 1;i <= fe.N.size();i++)
	for (size_t k = 1;k <= cmp;k++)
	  for (size_t l = 1;l <= cmp;l++)
	    // Diffusion
	    EV((i-1)*nf + k) -= scale*dUdX(l,k)*fe.dNdX(i,l)*fe.detJxW;
    }
  }
}

#endif
