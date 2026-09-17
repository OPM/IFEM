//==============================================================================
//!
//! \file PiolaOperators.C
//!
//! \date Apr 30 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete Piola-mapped operators.
//!
//==============================================================================

#include "PiolaOperators.h"
#include "FiniteElement.h"
#include "Tensor.h"
#include "Vec3.h"

#include <utility>


namespace {

void AdvectionConvInt (Matrix& C,
                       const FiniteElement& fe,
                       const Vec3& U, double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Matrix G(nsd, fe.dPdX.cols());
  for (size_t l = 1; l <= nsd; ++l)
  {
    for (size_t k = 1; k <= nsd; ++k)
      G.fillRow(k, fe.dPdX.getRow((l-1)*nsd+k).ptr());
    C.multiply(fe.P, G, true, false, l > 1, U[l-1] * scale * fe.detJxW);
  }
}


//! \brief Contracts the gradient of the mapped basis with a vector.
void GradientDotVector (Matrix& W,
                        const FiniteElement& fe,
                        const Vec3& U)
{
  const size_t nsd = fe.dNdX.cols();
  W.resize(nsd, fe.dPdX.cols(), true);
  for (size_t l = 1; l <= nsd; ++l)
    for (size_t k = 1; k <= nsd; ++k)
      for (size_t i = 1; i <= fe.dPdX.cols(); ++i)
        W(l,i) += U[k-1] * fe.dPdX((l-1)*nsd+k, i);
}


//! \brief The weights of the convective and the conservative term in a form.
std::pair<double,double> FormWeights (WeakOperators::ConvectionForm form)
{
  switch (form) {
    case WeakOperators::CONSERVATIVE:  return { 0.0,  -1.0};
    case WeakOperators::SKEWSYMMETRIC: return { 0.5,  -0.5};
    default:                           return { 1.0,   0.0};
  }
}

}


void PiolaOperators::Weak::Advection (Matrices& EM,
                                      const FiniteElement& fe,
                                      const Vec3& AC,
                                      const std::array<std::array<int,3>,3>& idx,
                                      double scale,
                                      WeakOperators::ConvectionForm cnvForm)
{
  Matrix A;
  AdvectionConvInt(A, fe, AC, scale);

  const auto [cw, tw] = FormWeights(cnvForm);
  Matrix C(A.rows(), A.cols());
  for (size_t i = 1; i <= C.rows(); ++i)
    for (size_t j = 1; j <= C.cols(); ++j)
      C(i,j) = cw*A(i,j) + tw*A(j,i);

  Copy(EM, fe, idx, C);
}


void PiolaOperators::Weak::Convection (Matrices& EM,
                                       const FiniteElement& fe,
                                       const Vec3& U,
                                       const Tensor& dUdX,
                                       const std::array<std::array<int,3>,3>& idx,
                                       double scale,
                                       WeakOperators::ConvectionForm form)
{
  Matrix A;
  AdvectionConvInt(A, fe, U, scale);

  Matrix dudx(dUdX.dim(), dUdX.dim());
  dudx = dUdX;
  Matrix dudxP;
  dudxP.multiply(dudx, fe.P);
  Matrix N;
  N.multiply(fe.P, dudxP, true, false, false, scale*fe.detJxW);

  Matrix W, M;
  GradientDotVector(W, fe, U);
  M.multiply(W, fe.P, true, false, false, scale*fe.detJxW);

  const auto [cw, tw] = FormWeights(form);
  Matrix C(A.rows(), A.cols());
  for (size_t i = 1; i <= C.rows(); ++i)
    for (size_t j = 1; j <= C.cols(); ++j)
      C(i,j) = cw*(A(i,j) + N(i,j)) + tw*(A(j,i) + M(i,j));

  Copy(EM, fe, idx, C);
}


void PiolaOperators::Weak::Gradient (Matrices& EM,
                                     const FiniteElement& fe,
                                     std::array<int,3>& idx,
                                     double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Vector divVel(fe.dPdX.cols());
  for (size_t i = 1; i <= fe.dPdX.cols(); ++i)
    for (size_t d = 1; d <= nsd; ++d)
      divVel(i) += fe.dPdX(1+(d-1)*(nsd+1),i);
  Matrix D;
  D.outer_product(divVel, fe.basis(nsd+1));
  for (size_t j = 1; j <= fe.basis(nsd+1).size(); ++j) {
    size_t ofs = 0;
    for (size_t b = 1; b <= nsd; ++b) {
      for (size_t i = 1; i <= fe.basis(b).size(); ++i)
        EM[idx[b-1]](i,j) -= D(i+ofs,j) * scale * fe.detJxW;
      ofs += fe.basis(b).size();
    }
  }
}


void PiolaOperators::Weak::ItgConstraint (std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          const std::array<int,3>& idx)
{
  const size_t nsd = fe.dNdX.cols();
  size_t ofs = 0;
  for (size_t b = 1; b <= nsd; ++b) {
    Matrix& EMb = EM[idx[b-1]];
    for (size_t i = 1; i <= fe.basis(b).size(); ++i)
      for (size_t k = 1; k <= nsd; ++k)
        EMb(i,k) += fe.P(k,ofs+i) * fe.detJxW;
    ofs += fe.basis(b).size();
  }
}


void PiolaOperators::Weak::Laplacian (Matrices& EM,
                                      const FiniteElement& fe,
                                      const std::array<std::array<int,3>,3>& idx,
                                      double scale, bool stress)
{
  Matrix A;
  A.multiply(fe.dPdX, fe.dPdX, true, false, false, scale*fe.detJxW);
  if (stress) {
    const size_t nsd = fe.dNdX.cols();
    Matrix dPdXT(fe.dPdX.rows(), fe.dPdX.cols());
    for (size_t d = 1; d <= nsd; ++d)
      for (size_t j = 1; j <= nsd; ++j) {
        const Vector row = fe.dPdX.getRow((j-1)*nsd + d);
        dPdXT.fillRow((d-1)*nsd + j, row.ptr());
      }
    A.multiply(fe.dPdX, dPdXT, true, false, true, scale*fe.detJxW);
  }
  Copy(EM, fe, idx, A);
}


bool PiolaOperators::Weak::MassCoeff (Matrices& EM,
                                      const Matrix& C,
                                      const FiniteElement& fe,
                                      const std::array<std::array<int,3>,3>& idx,
                                      double scale)
{
  Matrix CP, M;
  CP.multiply(C, fe.P);
  M.multiply(fe.P, CP, true, false, false, scale*fe.detJxW);
  Copy(EM, fe, idx, M);

  return true;
}


void PiolaOperators::Weak::Mass (Matrices& EM,
                                 const FiniteElement& fe,
                                 const std::array<std::array<int,3>,3>& idx,
                                 double scale)
{
  Matrix M;
  M.multiply(fe.P, fe.P, true, false, false, scale*fe.detJxW);
  Copy(EM, fe, idx, M);
}


void PiolaOperators::Weak::Source (Vectors& EV, const FiniteElement& fe,
                                   const Vec3& f, const std::array<int,3>& idx,
                                   double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Vector fV(f.ptr(), nsd);
  Vector Mv;
  fe.P.multiply(fV, Mv, scale * fe.detJxW, 0.0, true);
  Copy(EV, fe, idx, Mv);
}


void PiolaOperators::Residual::Convection (Vectors& EV, const FiniteElement& fe,
                                           const Vec3& U, const Tensor& dUdX,
                                           const Vec3& UC,
                                           const std::array<int,3>& idx, double scale,
                                           WeakOperators::ConvectionForm form)
{
  const size_t nsd = fe.grad(1).cols();
  const auto [cw, tw] = FormWeights(form);

  Vector r(fe.dPdX.cols());

  if (cw != 0.0) {
    Matrix C(nsd,1);
    C.fillColumn(1, (dUdX*UC).ptr());
    Matrix T;
    T.multiply(fe.P, C, true, false, false, -cw*scale*fe.detJxW);
    r.add(T, 1.0);
  }

  if (tw != 0.0) {
    Tensor UUC(nsd);
    for (size_t l = 1; l <= nsd; ++l)
      for (size_t k = 1; k <= nsd; ++k)
        UUC(k,l) = U[k-1] * UC[l-1];

    Vector diff;
    fe.dPdX.multiply(UUC, diff, true);
    r.add(diff, -tw*scale*fe.detJxW);
  }

  Copy(EV, fe, idx, r);
}


void PiolaOperators::Residual::Gradient (Vectors& EV,
                                         const FiniteElement& fe,
                                         const std::array<int,3>& idx, double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Vector divVel(fe.dPdX.cols());
  for (size_t i = 1; i <= fe.dPdX.cols(); ++i)
    for (size_t d = 1; d <= nsd; ++d)
      divVel(i) += fe.dPdX(1+(d-1)*(nsd+1),i) * scale * fe.detJxW;

  Copy(EV, fe, idx, divVel);
}


void PiolaOperators::Residual::Laplacian (Vectors& EV,
                                          const FiniteElement& fe,
                                          const Tensor& dUdX,
                                          const std::array<int,3>& idx,
                                          double scale, bool stress)
{
  Vector diff;
  fe.dPdX.multiply(dUdX, diff, true);
  if (stress) {
    Tensor dUdXT(dUdX);
    dUdXT.transpose();
    Vector diffT;
    fe.dPdX.multiply(dUdXT, diffT, true);
    diff += diffT;
  }
  diff *= -scale*fe.detJxW;
  Copy(EV, fe, idx, diff);
}


void PiolaOperators::Copy (Matrices& EM,
                           const FiniteElement& fe,
                           const std::array<std::array<int,3>,3>& idx,
                           const Matrix& A)
{
  const size_t nsd = fe.dNdX.cols();
  if (nsd < 1 || nsd > 3)
    return;
  size_t ofs = 1;
  for (size_t b = 1; b <= nsd; ++b) {
    size_t ofs2 = 1;
    for (size_t d = 1; d <= nsd; ++d) {
      if (!EM[idx[b-1][d-1]].empty())
        A.extractBlock(EM[idx[b-1][d-1]], ofs, ofs2, true);
      ofs2 += fe.basis(d).size();
    }
    ofs += fe.basis(b).size();
  }
}


void PiolaOperators::Copy (Vectors& EV,
                           const FiniteElement& fe,
                           const std::array<int,3>& idx,
                           const RealArray& V)
{
  size_t ofs = 0;
  const size_t nsd = fe.dNdX.cols();
  for (size_t b = 1; b <= nsd; ++b) {
    EV[idx[b-1]].add(V, 1.0, ofs);
    ofs += fe.basis(b).size();
  }
}
