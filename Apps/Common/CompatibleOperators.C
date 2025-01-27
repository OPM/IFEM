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
                                          const Vec3& AC,
                                          const std::array<std::array<int,3>,3>& idx,
                                          double scale,
                                          WeakOperators::ConvectionForm cnvForm)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    for (size_t k = 1; k <= nsd; ++k)
      if (cnvForm == WeakOperators::CONVECTIVE)
        EM[idx[n-1][n-1]].outer_product(fe.basis(n), fe.grad(n).getColumn(k),
                                        true,scale*AC[k-1]*fe.detJxW);
      else if (cnvForm == WeakOperators::CONSERVATIVE)
        EM[idx[n-1][n-1]].outer_product(fe.grad(n).getColumn(k),
                                        fe.basis(n), true, -scale*AC[k-1]*fe.detJxW);
      else if (cnvForm == WeakOperators::SKEWSYMMETRIC) {
        EM[idx[n-1][n-1]].outer_product(fe.basis(n), fe.grad(n).getColumn(k),
                                        true,0.5*scale*AC[k-1]*fe.detJxW);
        EM[idx[n-1][n-1]].outer_product(fe.grad(n).getColumn(k),
                                        fe.basis(n), true, -0.5*scale*AC[k-1]*fe.detJxW);
      }
}


void CompatibleOperators::Weak::Convection(std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const Vec3& U,
                                           const Tensor& dUdX,
                                           const std::array<std::array<int,3>,3>& idx,
                                           double scale,
                                           WeakOperators::ConvectionForm form)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t k = 1; k <= nsd; ++k)
    for (size_t l = 1; l <= nsd; ++l) {
      switch (form) {
        case WeakOperators::CONVECTIVE:
          EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.basis(l), true, scale*fe.detJxW*dUdX(k,l));
          if(k == l)
            for (size_t m = 1; m <= nsd; ++m)
              EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, scale*fe.detJxW*U[m-1]);
          break;
        case WeakOperators::CONSERVATIVE:
          EM[idx[k-1][l-1]].outer_product(fe.grad(k).getColumn(l), fe.basis(l), true, -U[k-1]*scale*fe.detJxW);
          if (k == l)
            for (size_t m = 1; m <= nsd; ++m)
              EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, -U[m-1]*scale*fe.detJxW);
          break;
        case WeakOperators::SKEWSYMMETRIC:
          EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.basis(l), true, 0.5*scale*fe.detJxW*dUdX(k,l));
          EM[idx[k-1][l-1]].outer_product(fe.grad(k).getColumn(l), fe.basis(l), true, -0.5*U[k-1]*scale*fe.detJxW);
          if (k == l)
            for (size_t m = 1; m <= nsd; ++m) {
              EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, 0.5*scale*fe.detJxW*U[m-1]);
              EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, -0.5*U[m-1]*scale*fe.detJxW);
            }
      }
    }
}


void CompatibleOperators::Weak::Gradient (std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          const std::array<int,3>& idx,
                                          double scale)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EM[idx[n-1]].outer_product(fe.grad(n).getColumn(n), fe.basis(nsd+1),
                               true, -scale*fe.detJxW);
}


void CompatibleOperators::Weak::Laplacian (std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const std::array<std::array<int,3>,3>& idx,
                                           double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EqualOrderOperators::Weak::Laplacian(EM[n], fe, scale, false, n);

  for (size_t m = 1; m <= nsd && stress; m++)
    for (size_t n = m; n <= nsd; n++) {
      const Vector mn = fe.grad(m).getColumn(n);
      const Vector nm = fe.grad(n).getColumn(m);
      EM[idx[m-1][n-1]].outer_product(mn, nm, true, scale*fe.detJxW);
      if (m != n && !EM[idx[n-1][m-1]].empty())
        EM[idx[n-1][m-1]].outer_product(nm, mn, true, scale*fe.detJxW);
    }
}


void CompatibleOperators::Weak::Mass (std::vector<Matrix>& EM,
                                      const FiniteElement& fe,
                                      const std::array<std::array<int,3>,3>& idx,
                                      double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Mass(EM[idx[k-1][k-1]], fe, scale, k);
}

void CompatibleOperators::Weak::Source (Vectors& EV,
                                        const FiniteElement& fe,
                                        const Vec3& f,
                                        const std::array<int, 3>& idx,
                                        double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[idx[k-1]], fe, scale*f[k-1], 1, k);
}


void CompatibleOperators::Weak::Source(Vectors& EV,
                                       const FiniteElement& fe,
                                       const std::array<int, 3>& idx,
                                       double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[idx[k-1]], fe, scale, 1, k);
}


void CompatibleOperators::Residual::Convection (Vectors& EV, const FiniteElement& fe,
                                                const Vec3& U, const Tensor& dUdX,
                                                const Vec3& UC, const std::array<int, 3>& idx,
                                                double scale, WeakOperators::ConvectionForm form)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t k = 1; k <= nsd; ++k) {
    double conv = 0.0;
    switch (form) {
      case WeakOperators::CONVECTIVE:
        for (size_t l = 1; l <= nsd; ++l)
          conv += -UC[l-1]*dUdX(k,l);
        EV[idx[k-1]].add(fe.basis(k), scale*conv*fe.detJxW);
        break;
     case WeakOperators::CONSERVATIVE:
        for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
          conv = 0.0;
          for (size_t l = 1; l <= nsd; ++l)
            conv += U[k-1]*UC[l-1]*fe.grad(k)(i,l);
          EV[idx[k-1]](i) += conv*scale*fe.detJxW;
        }
        break;
     case WeakOperators::SKEWSYMMETRIC:
        for (size_t l = 1; l <= nsd; ++l)
          conv += -UC[l-1]*dUdX(k,l);
        EV[idx[k-1]].add(fe.basis(k), 0.5*scale*conv*fe.detJxW);
        for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
          conv = 0.0;
          for (size_t l = 1; l <= nsd; ++l)
            conv += U[k-1]*UC[l-1]*fe.grad(k)(i,l);
          EV[idx[k-1]](i) += 0.5*scale*conv*fe.detJxW;
        }
        break;
    }
  }
}


void CompatibleOperators::Residual::Gradient (Vectors& EV, const FiniteElement& fe,
                                              const std::array<int, 3>& idx,
                                              double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EV[idx[k-1]].add(fe.grad(k).getColumn(k), scale*fe.detJxW);
}


void CompatibleOperators::Residual::Laplacian(Vectors& EV,
                                              const FiniteElement& fe,
                                              const Tensor& dUdX,
                                              const std::array<int, 3>& idx,
                                              double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  auto dUdXT = dUdX;
  dUdXT.transpose();
  for (size_t k = 1; k <= nsd; ++k) {
    Vector diff;
    fe.grad(k).multiply(Vector(dUdXT[k-1].ptr(), nsd), diff);
    if (stress)
      fe.grad(k).multiply(Vector(dUdX[k-1].ptr(), nsd), diff, false, 1);
    EV[idx[k-1]].add(diff, -scale*fe.detJxW);
  }
}
