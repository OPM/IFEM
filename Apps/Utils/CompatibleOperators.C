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
  const size_t nsd = fe.grad(1).cols();
  const double coef = scale*fe.detJxW;

  // The term is made of the piece which advects the trial function and the
  // piece which advects the test function
  const auto [cw, tw] = WeakOperators::FormConvectionWeights(form);

  // The advected function is of the component it belongs to, so unlike on an
  // equal order pair there is one such derivative for each of them
  std::vector<Vector> DU(nsd);
  for (size_t k = 1; k <= nsd; ++k) {
    DU[k-1].resize(fe.basis(k).size());
    for (size_t m = 1; m <= nsd; ++m)
      DU[k-1].add(fe.grad(k).getColumn(m), U[m-1]);
  }

  for (size_t k = 1; k <= nsd; ++k)
    for (size_t l = 1; l <= nsd; ++l) {
      Matrix& C = EM[idx[k-1][l-1]];
      if (cw != 0.0)
        C.outer_product(fe.basis(k), fe.basis(l), true, cw*coef*dUdX(k,l));
      if (tw != 0.0)
        C.outer_product(fe.grad(k).getColumn(l), fe.basis(l), true,
                        tw*coef*U[k-1]);
      if (k == l) {
        if (cw != 0.0)
          C.outer_product(fe.basis(k), DU[l-1], true, cw*coef);
        if (tw != 0.0)
          C.outer_product(DU[k-1], fe.basis(l), true, tw*coef);
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

void CompatibleOperators::Weak::ItgConstraint (std::vector<Matrix>& EM,
                                               const FiniteElement& fe,
                                              const std::array<int,3>& idx)
{
  const size_t nsd = fe.grad(1).cols();
  for (size_t i = 1; i <= nsd; ++i) {
    Matrix& EMi = EM[idx[i-1]];
    for (size_t j = 1; j <= fe.basis(i).size(); ++j)
      EMi(j,i) += fe.basis(i)(j) * fe.detJxW;
  }
}


void CompatibleOperators::Weak::Laplacian (std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const std::array<std::array<int,3>,3>& idx,
                                           double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EqualOrderOperators::Weak::Laplacian(EM[idx[n-1][n-1]], fe, scale, false, n);

  for (size_t m = 1; m <= nsd && stress; m++)
    for (size_t n = m; n <= nsd; n++) {
      const Vector mn = fe.grad(m).getColumn(n);
      const Vector nm = fe.grad(n).getColumn(m);
      EM[idx[m-1][n-1]].outer_product(mn, nm, true, scale*fe.detJxW);
      if (m != n && !EM[idx[n-1][m-1]].empty())
        EM[idx[n-1][m-1]].outer_product(nm, mn, true, scale*fe.detJxW);
    }
}


bool CompatibleOperators::Weak::MassCoeff (std::vector<Matrix>& EM,
                                           const Matrix& C,
                                           const FiniteElement& fe,
                                           const std::array<std::array<int,3>,3>& idx,
                                           double scale)
{
  const size_t nsd = fe.grad(1).cols();
  for (size_t k = 1; k <= nsd; ++k)
    for (size_t l = 1; l <= nsd; ++l)
      if (C(k,l) != 0.0 && !EM[idx[k-1][l-1]].empty())
        EM[idx[k-1][l-1]].outer_product(fe.basis(k), fe.basis(l), true,
                                        scale*fe.detJxW*C(k,l));

  return true;
}


void CompatibleOperators::Weak::LaplacianCoeff (std::vector<Matrix>& EM,
                                                const Matrix& K,
                                                const FiniteElement& fe,
                                                const std::array<std::array<int,3>,3>& idx,
                                                double scale)
{
  const size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EqualOrderOperators::Weak::LaplacianCoeff(EM[idx[n-1][n-1]], K, fe, scale, n);
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


void CompatibleOperators::Weak::Source (Vectors& EV,
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
  const size_t nsd = fe.grad(1).cols();
  const double coef = scale*fe.detJxW;

  // The term is made of the piece which advects the velocity and the piece
  // which advects the test function
  const auto [cw, tw] = WeakOperators::FormConvectionWeights(form);

  for (size_t k = 1; k <= nsd; ++k) {
    const Vector& N = fe.basis(k);
    const Matrix& D = fe.grad(k);
    for (size_t i = 1; i <= N.size(); ++i)
      for (size_t l = 1; l <= nsd; ++l) {
        double conv = 0.0;
        if (cw != 0.0) conv -= cw*UC[l-1]*dUdX(k,l)*N(i);
        if (tw != 0.0) conv -= tw*U[k-1]*UC[l-1]*D(i,l);
        EV[idx[k-1]](i) += coef*conv;
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
