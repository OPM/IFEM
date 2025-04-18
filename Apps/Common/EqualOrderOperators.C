//==============================================================================
//!
//! \file EqualOrderOperators.C
//!
//! \date Jul 22 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various discrete equal-ordered operators.
//!
//==============================================================================

#include "EqualOrderOperators.h"
#include "FiniteElement.h"
#include "Vec3.h"
#include "Vec3Oper.h"

namespace {

//! \brief Helper for adding an element matrix to one component.
//! \param[out] EM The element matrix to add to.
//! \param[in] A The scalar element matrix to add.
//! \param[in] cmp1 First components to add matrix to
//! \param[in] cmp2 Second irst components to add matrix to
//! \param[in] nf Number of components in total matrix.
void addComponent (Matrix& EM, const Matrix& A,
                   size_t cmp1, size_t cmp2, size_t nf)
{
  for (size_t i = 1; i <= A.rows(); ++i)
    for (size_t j = 1; j <= A.cols(); ++j)
      EM(nf*(i-1)+cmp1,nf*(j-1)+cmp2) += A(i, j);
}


//! \brief Helper for adding an element matrix to several components.
//! \param[out] EM The element matrix to add to.
//! \param[in] A The scalar element matrix to add.
//! \param[in] cmp Number of components to add matrix to
//! \param[in] nf Number of components in total matrix.
//! \param[in] scmp Index of first component to add matrix to.
void addComponents (Matrix& EM, const Matrix& A,
                    size_t cmp, size_t nf, size_t scmp)
{
  if (cmp == 1 && nf == 1)
    EM += A;
  else
    for (size_t k = 1; k <= cmp; ++k)
      addComponent(EM, A, k+scmp, k+scmp, nf);
}


//! \brief Helper applying a divergence (1) or a gradient (2) operation
template<int Operation>
void DivGrad (Matrix& EM, const FiniteElement& fe,
              double scale, int basis, int tbasis)
{
  size_t nsd = fe.grad(basis).cols();
  for (size_t i = 1; i <= fe.basis(tbasis).size();++i)
    for (size_t j = 1; j <= fe.basis(basis).size();++j)
      for (size_t k = 1; k <= nsd; ++k) {
        double div = fe.basis(basis)(j)*fe.grad(tbasis)(i,k)*fe.detJxW;
        if (Operation == 2)
          EM((i-1)*nsd+k,j) += -scale*div;
        if (Operation == 1)
          EM(j, (i-1)*nsd+k) += scale*div;
      }
}

}


void EqualOrderOperators::Weak::Advection (Matrix& EM, const FiniteElement& fe,
                                           const Vec3& AC,
                                           double scale,
                                           WeakOperators::ConvectionForm form,
                                           int basis)
{
  Matrix C(fe.basis(basis).size(), fe.basis(basis).size());
  size_t ncmp = EM.rows() / C.rows();

  // Sum convection for each direction
  for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
    if (form == WeakOperators::CONVECTIVE)
      C.outer_product(fe.basis(basis),
                      fe.grad(basis).getColumn(k), true,
                      scale*AC[k-1]*fe.detJxW);
    else if (form == WeakOperators::CONSERVATIVE)
      C.outer_product(fe.grad(basis).getColumn(k),
                      fe.basis(basis), true,
                      -scale*AC[k-1]*fe.detJxW);
    else if (form == WeakOperators::SKEWSYMMETRIC) {
      C.outer_product(fe.basis(basis),
                      fe.grad(basis).getColumn(k), true,
                      0.5*scale*AC[k-1]*fe.detJxW);
      C.outer_product(fe.grad(basis).getColumn(k),
                      fe.basis(basis), true,
                      -0.5*scale*AC[k-1]*fe.detJxW);
    }

  addComponents(EM, C, ncmp, ncmp, 0);
}


void EqualOrderOperators::Weak::Convection (Matrix& EM, const FiniteElement& fe,
                                            const Vec3& U, const Tensor& dUdX,
                                            double scale,
                                            WeakOperators::ConvectionForm form,
                                            int basis)
{
  size_t cmp = EM.rows() / fe.basis(basis).size();
  double coef = scale*fe.detJxW;

  const Vector& N = fe.basis(basis);
  const Matrix& D = fe.grad(basis);

  Matrix B;
  if (form != WeakOperators::CONSERVATIVE)
    B.outer_product(N, N);
  for (size_t k = 1; k <= cmp; ++k)
    for (size_t l = 1; l <= cmp; ++l) {
      Matrix C(N.size(), N.size());
      switch (form) {
        case WeakOperators::CONVECTIVE:
          C.add(B, dUdX(k,l));
          if (k == l)
            for (size_t m = 1; m <= cmp; ++m)
              C.outer_product(N, D.getColumn(m), true, U[m-1]);
          break;
        case WeakOperators::CONSERVATIVE:
          C.outer_product(D.getColumn(l), N, false, -U[k-1]);
          if (k == l)
            for (size_t m = 1; m <= cmp; ++m)
              C.outer_product(D.getColumn(m), N, true, -U[m-1]);
          break;
        case WeakOperators::SKEWSYMMETRIC:
          C.add(B, 0.5*dUdX(k,l));
          C.outer_product(D.getColumn(l), N, true, -0.5*U[k-1]);
          if (k == l)
            for (size_t m = 1; m <= cmp; ++m) {
              C.outer_product(N, D.getColumn(m), true, 0.5*U[m-1]);
              C.outer_product(D.getColumn(m), N, true, -0.5*U[m-1]);
            }
          break;
      }
      for (size_t i = 1; i <= N.size(); i++)
        for (size_t j = 1; j <= N.size(); j++)
          EM((i-1)*cmp+k,(j-1)*cmp+l) += coef*C(i,j);
    }
}


void EqualOrderOperators::Weak::Divergence (Matrix& EM, const FiniteElement& fe,
                                            double scale, int basis, int tbasis)
{
  DivGrad<1>(EM,fe,scale,basis,tbasis);
}


void EqualOrderOperators::Weak::Gradient (Matrix& EM, const FiniteElement& fe,
                                          double scale, int basis, int tbasis)
{
  DivGrad<2>(EM,fe,scale,basis,tbasis);
}


void EqualOrderOperators::Weak::Divergence (Vector& EV,
                                            const FiniteElement& fe,
                                            const Vec3& D,
                                            double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  fe.grad(basis).multiply(Vector(D.ptr(),nsd),EV,scale*fe.detJxW,1.0);
}


void EqualOrderOperators::Weak::Laplacian (Matrix& EM, const FiniteElement& fe,
                                           double scale, bool stress, int basis)
{
  size_t cmp = EM.rows() / fe.basis(basis).size();
  Matrix A;
  A.multiply(fe.grad(basis),fe.grad(basis),false,true,false,scale*fe.detJxW);
  addComponents(EM, A, cmp, cmp, 0);
  if (stress)
    Stress(EM, fe, scale, basis);
}


void EqualOrderOperators::Weak::Stress (Matrix& EM,
                                        const FiniteElement& fe,
                                        double scale, int basis)
{
  Matrix A;
  size_t cmp = EM.rows() / fe.basis(basis).size();
  for (size_t k = 1; k <= cmp; k++)
    for (size_t l = 1; l <= cmp; l++) {
      A.outer_product(fe.grad(basis).getColumn(l),
                      fe.grad(basis).getColumn(k),
                      false, scale * fe.detJxW);
      addComponent(EM, A, k, l, cmp);
    }
}


void EqualOrderOperators::Weak::LaplacianCoeff (Matrix& EM, const Matrix& K,
                                                const FiniteElement& fe,
                                                double scale, int basis)
{
  Matrix KB;
  KB.multiply(K,fe.grad(basis),false,true).multiply(scale*fe.detJxW);
  EM.multiply(fe.grad(basis),KB,false,false,true);
}


void EqualOrderOperators::Weak::Mass (Matrix& EM, const FiniteElement& fe,
                                      double scale, int basis)
{
  size_t ncmp = EM.rows()/fe.basis(basis).size();
  Matrix A;
  A.outer_product(fe.basis(basis),fe.basis(basis),false);
  A *= scale*fe.detJxW;
  addComponents(EM, A, ncmp, ncmp, 0);
}


void EqualOrderOperators::Weak::ItgConstraint (Matrix& EM, const FiniteElement& fe,
                                               double scale, int basis)
{
  const size_t ncmp = EM.rows() / fe.basis(basis).size();
  Matrix A(fe.basis(basis).size(), 1);
  A = fe.basis(basis);
  A *= scale*fe.detJxW;
  addComponents(EM, A, ncmp, ncmp, 0);
}


void EqualOrderOperators::Weak::Source (Vector& EV,
                                        const FiniteElement& fe,
                                        double scale, int cmp, int basis)
{
  size_t ncmp = EV.size() / fe.basis(basis).size();
  if (cmp == 1 && ncmp == 1)
    EV.add(fe.basis(basis), scale*fe.detJxW);
  else {
    for (size_t k  = (cmp == 0 ? 1 : cmp);
                k <= (cmp == 0 ? ncmp : cmp); ++k)
      EV.add(fe.basis(basis), scale*fe.detJxW, 0, 1, k-1, ncmp);
  }
}


void EqualOrderOperators::Weak::Source (Vector& EV, const FiniteElement& fe,
                                        const Vec3& f, double scale, int basis)
{
  size_t cmp = EV.size() / fe.basis(basis).size();
  for (size_t k = 1; k <= cmp; ++k)
    EV.add(fe.basis(basis), scale*f[k-1]*fe.detJxW, 0, 1, k-1, cmp);
}


void EqualOrderOperators::Residual::Advection (Vector& EV,
                                               const FiniteElement& fe,
                                               const Vec3& AC, const Tensor& g,
                                               double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  for (size_t k = 1; k <= nsd; ++k)
    EV.add(fe.basis(basis), (g[k-1]*AC)*scale*fe.detJxW, 0, 1, k-1, nsd);
}


void EqualOrderOperators::Residual::Convection (Vector& EV,
                                                const FiniteElement& fe,
                                                const Vec3& U,
                                                const Tensor& dUdX,
                                                const Vec3& UC, double scale,
                                                WeakOperators::ConvectionForm form,
                                                int basis)
{
  size_t cmp = EV.size() / fe.basis(basis).size();
  double coef = scale * fe.detJxW;
  double conv = 0.0;
  for (size_t i = 1;i <= fe.basis(basis).size();i++)
    for (size_t k = 1;k <= cmp;k++)
      for (size_t l = 1;l <= cmp;l++) {
        switch (form) {
          case WeakOperators::CONVECTIVE:
            conv = -UC[l-1]*dUdX(k,l)*fe.basis(basis)(i);
            break;
          case WeakOperators::CONSERVATIVE:
            conv = U[k-1]*UC[l-1]*fe.grad(basis)(i,l);
            break;
          case WeakOperators::SKEWSYMMETRIC:
            conv = U[k-1]*UC[l-1]*fe.grad(basis)(i,l)
                 - UC[l-1]*dUdX(k,l)*fe.basis(basis)(i);
            conv *= 0.5;
            break;
          default:
            std::cerr << "EqualOrderOperators::Residual::Convection: "
                      << "Unknown form " << form << std::endl;
        }
        EV((i-1)*cmp+k) += coef*conv;
      }
}


void EqualOrderOperators::Residual::Divergence (Vector& EV,
                                                const FiniteElement& fe,
                                                const Tensor& dUdX,
                                                double scale, size_t basis)
{
  EV.add(fe.basis(basis),-scale*dUdX.trace()*fe.detJxW);
}


void EqualOrderOperators::Residual::Gradient (Vector& EV,
                                              const FiniteElement& fe,
                                              double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  for (size_t k = 1; k <= nsd; ++k)
    EV.add(fe.grad(basis).getColumn(k), scale*fe.detJxW, 0, 1, k-1, nsd);
}




void EqualOrderOperators::Residual::Laplacian (Vector& EV,
                                               const FiniteElement& fe,
                                               const Vec3& dUdX,
                                               double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  fe.grad(basis).multiply(Vector(dUdX.ptr(),nsd),EV,-scale*fe.detJxW,1.0);
}

void EqualOrderOperators::Residual::Laplacian (Vector& EV,
                                               const FiniteElement& fe,
                                               const Tensor& dUdX,
                                               double scale,
                                               bool stress, int basis)
{
  size_t nsd = fe.grad(1).cols();
  auto dUdXT = dUdX;
  dUdXT.transpose();
  for (size_t k = 1; k <= nsd; ++k) {
    Vector diff;
    fe.grad(basis).multiply(Vector(dUdXT[k-1].ptr(), nsd), diff);
    if (stress)
      fe.grad(basis).multiply(Vector(dUdX[k-1].ptr(), nsd), diff, false, 1);
    EV.add(diff, -scale*fe.detJxW, 0, 1, k-1, nsd);
  }
}
