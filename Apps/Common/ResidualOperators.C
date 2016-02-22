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
  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Matrix& dUdX, double scale, size_t cmp,
                  size_t nf, size_t basis)
  {
    for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
      double div=0.0;
      for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
        div += dUdX(k,k);
      EV((i-1)*nf+cmp) += scale*div*fe.basis(basis)(i)*fe.detJxW;
    }
  }
}
