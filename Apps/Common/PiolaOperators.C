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
#include "Vec3.h"


namespace {

//! \brief Block indices for velocity blocks.
static constexpr int vidx[3][3] = {{ 1,  7,  8},
                                     {12,  2, 13},
                                     {17, 18,  3}};


void AdvectionConvInt (Matrix& C,
                       const FiniteElement& fe,
                       const Vec3& U, double scale)
{
  Matrix G1(2, fe.dPdX.cols()), G2(2, fe.dPdX.cols());
  G1.fillRow(1, fe.dPdX.getRow(1).data());
  G1.fillRow(2, fe.dPdX.getRow(3).data());
  G2.fillRow(1, fe.dPdX.getRow(2).data());
  G2.fillRow(2, fe.dPdX.getRow(4).data());
  C.multiply(fe.P, G1, true, false, false, U[0] * scale * fe.detJxW);
  C.multiply(fe.P, G2, true, false, true,  U[1] * scale * fe.detJxW);
}

}


void PiolaOperators::Weak::Advection (Matrices& EM,
                                      const FiniteElement& fe,
                                      const Vec3& AC, double scale,
                                      WeakOperators::ConvectionForm cnvForm)
{
  if (cnvForm != WeakOperators::CONVECTIVE) {
    std::cerr << "Only convective advection operator implemented with piola" << std::endl;
    exit(1);
  }

  Matrix C;
  AdvectionConvInt(C, fe, AC, scale);
  Copy(EM, fe, C);
}


void PiolaOperators::Weak::Convection (Matrices& EM,
                                       const FiniteElement& fe,
                                       const Vec3& U,
                                       const Tensor& dUdX, double scale,
                                       WeakOperators::ConvectionForm form)
{
  if (form != WeakOperators::CONVECTIVE) {
    std::cerr << "Only convective convection operator implemented with piola" << std::endl;
    exit(1);
  }

  Matrix C;
  Matrix dudx(dUdX.dim(), dUdX.dim());
  dudx = dUdX;

  Matrix C1;
  C1.multiply(dudx, fe.P);
  AdvectionConvInt(C, fe, U, scale);
  C.multiply(fe.P, C1, true, false, true, scale*fe.detJxW);
  Copy(EM, fe, C);
}


void PiolaOperators::Weak::Gradient (Matrices& EM,
                                     const FiniteElement& fe,
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
        EM[9 + 5*(b-1)](i,j) -= D(i+ofs,j) * scale * fe.detJxW;
      ofs += fe.basis(b).size();
    }
  }
}


void PiolaOperators::Weak::Laplacian (Matrices& EM,
                                      const FiniteElement& fe,
                                      double scale, bool stress)
{
  if (stress) {
    std::cerr << "Stress laplacian operator not implemented with piola" << std::endl;
    exit(1);
  }
  Matrix A;
  A.multiply(fe.dPdX, fe.dPdX, true, false, false, scale*fe.detJxW);
  Copy(EM, fe, A);
}


void PiolaOperators::Weak::Mass (Matrices& EM,
                                 const FiniteElement& fe, double scale)
{
  Matrix M;
  M.multiply(fe.P, fe.P, true, false, false, scale*fe.detJxW);
  Copy(EM, fe, M);
}


void PiolaOperators::Weak::Source (Vectors& EV, const FiniteElement& fe,
                                   const Vec3& f, double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Vector fV(f.ptr(), nsd);
  Vector Mv;
  fe.P.multiply(fV, Mv, scale * fe.detJxW, 0.0, true);
  Copy(EV, fe, Mv);
}


void PiolaOperators::Residual::Convection (Vectors& EV, const FiniteElement& fe,
                                           const Vec3& U, const Tensor& dUdX,
                                           const Vec3& UC, double scale,
                                           WeakOperators::ConvectionForm form)
{
  if (form != WeakOperators::CONVECTIVE) {
    std::cerr << "Only convective convection implemented with piola" << std::endl;
     exit(1);
  }

  size_t nsd = fe.grad(1).cols();
  Matrix C(nsd,1);
  C.fillColumn(1, (dUdX*UC).ptr());
  Matrix T;
  T.multiply(fe.P, C, true, false, false, -scale*fe.detJxW);
  Copy(EV, fe, T);
}


void PiolaOperators::Residual::Gradient (Vectors& EV,
                                         const FiniteElement& fe,
                                         double scale)
{
  const size_t nsd = fe.dNdX.cols();
  Vector divVel(fe.dPdX.cols());
  for (size_t i = 1; i <= fe.dPdX.cols(); ++i)
    for (size_t d = 1; d <= nsd; ++d)
      divVel(i) += fe.dPdX(1+(d-1)*(nsd+1),i) * scale * fe.detJxW;

  Copy(EV, fe, divVel);
}


void PiolaOperators::Residual::Laplacian (Vectors& EV,
                                          const FiniteElement& fe,
                                          const Tensor& dUdX,
                                          double scale, bool stress)
{
  if (stress) {
    std::cerr << "Stress laplacian operator not implemented with piola" << std::endl;
    exit(1);
  }
  Vector diff;
  fe.dPdX.multiply(dUdX, diff, true);
  diff *= -scale*fe.detJxW;
  Copy(EV, fe, diff);
}


void PiolaOperators::Copy (Matrices& EM,
                           const FiniteElement& fe,
                           const Matrix& A)
{
  const size_t nsd = fe.dNdX.cols();
  if (nsd < 1 || nsd > 3)
      return;
  size_t ofs = 1;
  for (size_t b = 1; b <= nsd; ++b) {
    size_t ofs2 = ofs;
    for (size_t d = b; d <= nsd; ++d) {
      if (!EM[vidx[b-1][d-1]].empty())
        A.extractBlock(EM[vidx[b-1][d-1]], ofs, ofs2, true);
      ofs2 += fe.basis(d).size();
    }
    ofs += fe.basis(b).size();
  }
}


void PiolaOperators::Copy (Vectors& EV,
                           const FiniteElement& fe,
                           const Vector& V)
{
  size_t ofs = 0;
  const size_t nsd = fe.dNdX.cols();
  for (size_t b = 1; b <= nsd; ++b) {
    EV[b].add(V, 1.0, ofs);
    ofs += fe.basis(b).size();
  }
}
