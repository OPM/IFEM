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


namespace {
  //! \brief Block indices for velocity blocks.
  static constexpr int vidx[3][3] = {{ 1,  7,  8},
                                     {12,  2, 13},
                                     {17, 18,  3}};
}


void CompatibleOperators::Weak::Advection(std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          const Vec3& AC, double scale,
                                          WeakOperators::ConvectionForm cnvForm)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    for (size_t k = 1; k <= nsd; ++k)
      if (cnvForm == WeakOperators::CONVECTIVE)
        EM[n].outer_product(fe.basis(n), fe.grad(n).getColumn(k),
                            true,scale*AC[k-1]*fe.detJxW);
      else if (cnvForm == WeakOperators::CONSERVATIVE)
        EM[n].outer_product(fe.grad(n).getColumn(k),
                            fe.basis(n), true, -scale*AC[k-1]*fe.detJxW);
      else if (cnvForm == WeakOperators::SKEWSYMMETRIC) {
        EM[n].outer_product(fe.basis(n), fe.grad(n).getColumn(k),
                            true,0.5*scale*AC[k-1]*fe.detJxW);
        EM[n].outer_product(fe.grad(n).getColumn(k),
                            fe.basis(n), true, -0.5*scale*AC[k-1]*fe.detJxW);
      }
}


void CompatibleOperators::Weak::Convection(std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const Vec3& U,
                                           const Tensor& dUdX,
                                           double scale,
                                           WeakOperators::ConvectionForm form)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t k = 1; k <= nsd; ++k)
    for (size_t l = 1; l <= nsd; ++l) {
      switch (form) {
        case WeakOperators::CONVECTIVE:
          EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.basis(l), true, scale*fe.detJxW*dUdX(k,l));
          if(k == l)
            for (size_t m = 1; m <= nsd; ++m)
              EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, scale*fe.detJxW*U[m-1]);
          break;
        case WeakOperators::CONSERVATIVE:
          EM[vidx[k-1][l-1]].outer_product(fe.grad(k).getColumn(l), fe.basis(l), true, -U[k-1]*scale*fe.detJxW);
          if (k == l)
            for (size_t m = 1; m <= nsd; ++m)
              EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, -U[m-1]*scale*fe.detJxW);
          break;
        case WeakOperators::SKEWSYMMETRIC:
          EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.basis(l), true, 0.5*scale*fe.detJxW*dUdX(k,l));
          EM[vidx[k-1][l-1]].outer_product(fe.grad(k).getColumn(l), fe.basis(l), true, -0.5*U[k-1]*scale*fe.detJxW);
          if (k == l)
            for (size_t m = 1; m <= nsd; ++m) {
              EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, 0.5*scale*fe.detJxW*U[m-1]);
              EM[vidx[k-1][l-1]].outer_product(fe.basis(k), fe.grad(l).getColumn(m), true, -0.5*U[m-1]*scale*fe.detJxW);
            }
      }
    }
}


void CompatibleOperators::Weak::Gradient(std::vector<Matrix>& EM,
                                        const FiniteElement& fe,
                                        double scale)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EM[9+5*(n-1)].outer_product(fe.grad(n).getColumn(n), fe.basis(nsd+1),
                                true, -scale*fe.detJxW);
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
      const Vector mn = fe.grad(m).getColumn(n);
      const Vector nm = fe.grad(n).getColumn(m);
      EM[vidx[m-1][n-1]].outer_product(mn, nm, true, scale*fe.detJxW);
      if (m != n && !EM[vidx[n-1][m-1]].empty())
        EM[vidx[n-1][m-1]].outer_product(nm, mn, true, scale*fe.detJxW);
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
  for (size_t k = 1; k <= nsd; ++k) {
    double conv = 0.0;
    switch (form) {
      case WeakOperators::CONVECTIVE:
        for (size_t l = 1; l <= nsd; ++l)
          conv += -UC[l-1]*dUdX(k,l);
        EV[k].add(fe.basis(k), scale*conv*fe.detJxW);
        break;
     case WeakOperators::CONSERVATIVE:
        for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
          conv = 0.0;
          for (size_t l = 1; l <= nsd; ++l)
            conv += U[k-1]*UC[l-1]*fe.grad(k)(i,l);
          EV[k](i) += conv*scale*fe.detJxW;
        }
        break;
     case WeakOperators::SKEWSYMMETRIC:
        for (size_t l = 1; l <= nsd; ++l)
          conv += -UC[l-1]*dUdX(k,l);
        EV[k].add(fe.basis(k), 0.5*scale*conv*fe.detJxW);
        for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
          conv = 0.0;
          for (size_t l = 1; l <= nsd; ++l)
            conv += U[k-1]*UC[l-1]*fe.grad(k)(i,l);
          EV[k](i) += 0.5*scale*conv*fe.detJxW;
        }
        break;
    }
  }
}


void CompatibleOperators::Residual::Gradient (Vectors& EV, const FiniteElement& fe,
                                              double scale)
{
  for (size_t n = 1; n <= fe.grad(n).cols(); ++n)
    EV[n].add(fe.grad(n).getColumn(n), scale*fe.detJxW);
}


void CompatibleOperators::Residual::Laplacian(Vectors& EV,
                                              const FiniteElement& fe,
                                              const Tensor& dUdX,
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
    EV[k].add(diff, -scale*fe.detJxW);
  }
}
