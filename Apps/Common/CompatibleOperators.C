//==============================================================================
//!
//! \file CompatibleOperators.C
//!
//! \date Oct 9 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete div-compatible operators.
//!
//==============================================================================

#include "CompatibleOperators.h"
#include "EqualOrderOperators.h"
#include "FiniteElement.h"
#include "Vec3.h"


void CompatibleOperators::Weak::Advection(std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          const Vec3& AC, double scale)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    for (size_t k = 1; k <= nsd; ++k)
        EM[n].outer_product(fe.basis(n), fe.grad(n).getColumn(k),
                            true,scale*AC[k]*fe.detJxW);
}


void CompatibleOperators::Weak::Convection(std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const Vec3& U,
                                           const Tensor& dUdX,
                                           double scale,
                                           WeakOperators::ConvectionForm form)
{
  static const double vidx[3][3] = {{ 1,  7,  8},
                                    {12,  2, 13},
                                    {17, 18,  3}};
  size_t nsd = fe.grad(1).cols();
  for (size_t m=1; m <= nsd; ++m)
    for (size_t i=1; i <= fe.basis(m).size(); ++i)
      for (size_t n = 1; n <= nsd; ++n)
        for (size_t j=1; j <= fe.basis(n).size(); ++j) {
          double conv = fe.basis(n)(j) * dUdX(m,n);
          if (m == n)
            for (size_t k = 1; k <= nsd; ++k)
              conv += U[k-1] * fe.grad(n)(j,k);
          EM[vidx[m-1][n-1]](i,j) += scale * conv * fe.basis(m)(i) * fe.detJxW;
        }
}


void CompatibleOperators::Weak::Gradient(std::vector<Matrix>& EM,
                                        const FiniteElement& fe,
                                        double scale)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EM[10+5*(n-1)].outer_product(fe.grad(n).getColumn(n), fe.basis(nsd+1),
                                 true, -scale*fe.detJxW);
}


void CompatibleOperators::Weak::Gradient(Vectors& EV, const FiniteElement& fe,
                                         double scale)
{
  for (size_t n = 1; n <= fe.grad(n).cols(); ++n)
    EV[n].add(fe.grad(n).getColumn(n), scale*fe.detJxW);
}


void CompatibleOperators::Weak::Laplacian(std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EqualOrderOperators::Weak::Laplacian(EM[n], fe, scale, false, n);

  for (size_t m = 1; m <= nsd && stress; m++)
    for (size_t n = m; n <= nsd; n++) {
      int idx = m == n ? m : (m == 1 ? 6+n-m : 13);
      EM[idx].outer_product(fe.grad(m).getColumn(n),
                            fe.grad(n).getColumn(m), true, scale*fe.detJxW);
    }
}


void CompatibleOperators::Weak::Mass(std::vector<Matrix>& EM,
                                     const FiniteElement& fe, double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Mass(EM[k], fe, scale, k);
}

void CompatibleOperators::Weak::Source(Vectors& EV,
                                       const FiniteElement& fe,
                                       const Vec3& f, double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[k], fe, scale*f[k-1], 1, k);
}


void CompatibleOperators::Weak::Source(Vectors& EV,
                                       const FiniteElement& fe,
                                       double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[k], fe, scale, 1, k);
}


void CompatibleOperators::Residual::Convection(Vectors& EV, const FiniteElement& fe,
                                               const Vec3& U, const Tensor& dUdX,
                                               const Vec3& UC, double scale,
                                               WeakOperators::ConvectionForm form)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n=1; n <= nsd; ++n) {
    double conv = 0.0;
    for (size_t k = 1; k<= nsd; ++k)
      conv += U[k-1]*dUdX(n, k);
    EV[n].add(fe.basis(n), -scale*conv*fe.detJxW);
  }
}


void CompatibleOperators::Residual::Laplacian(Vectors& EV,
                                              const FiniteElement& fe,
                                              const Tensor& dUdX,
                                              double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t k = 1; k <= nsd; ++k)
    for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
      double diff = 0.0;
      for (size_t m = 1; m <= nsd; ++m)
        diff += fe.grad(k)(i,m)*dUdX(k,m);
      EV[k](i) += scale*diff*fe.detJxW;
    }
}
