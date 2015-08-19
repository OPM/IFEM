//==============================================================================
//!
//! \file ResidualOperators.C
//!
//! \date Aug 19 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete residual operators.
//!
//==============================================================================


#include "ResidualOperators.h"
#include "FiniteElement.h"
#include "Vec3.h"


namespace ResidualOperators {
  void Convection(Vector& EV, const FiniteElement& fe,
                  const Vec3& U, const Matrix& dUdX, const Vec3& UC,
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
      for (size_t i = 1;i <= fe.N.size();i++) {
        for (size_t k = 1;k <= cmp;k++) {
          // Convection
          double conv = 0.0;
          for (size_t l = 1;l <= cmp;l++)
            conv += UC[l-1]*dUdX(k,l);
          conv *= scale*fe.detJxW;

          EV((i-1)*nf + k) -= conv*fe.N(i);
        }
      }
    }
  }


  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Matrix& dUdX, double scale, size_t cmp, size_t nf)
  {
    for (size_t i = 1; i <= fe.N.size(); ++i) {
      double div=0.0;
      for (size_t k = 1; k <= fe.dNdX.cols(); ++k)
        div += dUdX(k,k);
      EV((i-1)*nf+cmp) += scale*div*fe.N(i)*fe.detJxW;
    }
  }


  void Laplacian(Vector& EV, const FiniteElement& fe,
                 const Matrix& dUdX, double scale, size_t cmp, size_t nf,
                 bool stress, size_t scmp)
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
