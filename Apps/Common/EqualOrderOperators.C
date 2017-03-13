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

//! \brief Helper for adding an element matrix to several components.
//! \param[out] EM The element matrix to add to.
//! \param[in] A The scalar element matrix to add.
//! \param[in] cmp Number of components to add matrix to
//! \param[in] nf Number of components in total matrix.
//! \param[in] scmp Index of first component to add matrix to.
static void addComponents(Matrix& EM, const Matrix& A,
                          size_t cmp, size_t nf, size_t scmp)
{
  if (cmp == 1 && nf == 1)
    EM += A;
  else
    for (size_t i = 1; i <= A.rows(); ++i)
      for (size_t j = 1; j <= A.cols(); ++j)
        for (size_t k = 1; k <= cmp; ++k)
          EM(nf*(i-1)+k+scmp,nf*(j-1)+k+scmp) += A(i, j);
}


//! \brief Helper applying a divergence (1) or a gradient (2) operation
template<int Operation>
static void DivGrad(Matrix& EM, const FiniteElement& fe,
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


void EqualOrderOperators::Weak::Advection(Matrix& EM, const FiniteElement& fe,
                                          const Vec3& AC,
                                          double scale, int basis)
{
  Matrix C(fe.basis(basis).size(), fe.basis(basis).size());
  size_t ncmp = EM.rows() / C.rows();
  for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
    for (size_t j = 1; j <= fe.basis(basis).size(); ++j) {
      // Sum convection for each direction
      for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
        C(i,j) += AC[k-1]*fe.grad(basis)(j,k);

      C(i,j) *= scale*fe.basis(basis)(i)*fe.detJxW;
    }
  }
  addComponents(EM, C, ncmp, ncmp, 0);
}


void EqualOrderOperators::Weak::Convection(Matrix& EM, const FiniteElement& fe,
                                           const Vec3& U, const Tensor& dUdX,
                                           double scale, WeakOperators::ConvectionForm form, int basis)
{
  size_t cmp = EM.rows() / fe.basis(basis).size();
  double coef = scale*fe.detJxW;
  double conv = 0.0;
  for (size_t i = 1;i <= fe.basis(basis).size();i++)
    for (size_t j = 1;j <= fe.basis(basis).size();j++)
      for (size_t k = 1;k <= cmp;k++)
        for (size_t l = 1;l <= cmp;l++) {
          switch (form) {
            case WeakOperators::CONVECTIVE:
              conv = fe.basis(basis)(i)*dUdX(k,l)*fe.basis(basis)(j);
              if (k==l)
                for (size_t m = 1;m <= cmp;m++)
                  conv += U[m-1]*fe.grad(basis)(j,m)*fe.basis(basis)(i);
              break;
            case WeakOperators::CONSERVATIVE:
              conv = -U[k-1]*fe.basis(basis)(j)*fe.grad(basis)(i,l);
              if (k==l)
                for (size_t m = 1;m <= cmp;m++)
                  conv -= fe.basis(basis)(j)*U[m-1]*fe.grad(basis)(i,m);
              break;
            case WeakOperators::SKEWSYMMETRIC:
              conv = fe.basis(basis)(i)*dUdX(k,l)*fe.basis(basis)(j)
                   - U[k-1]*fe.basis(basis)(j)*fe.grad(basis)(i,l);
              if (k==l)
                for (size_t m = 1;m <= cmp;m++)
                  conv += U[m-1]*fe.grad(basis)(j,m)*fe.basis(basis)(i)
                        - fe.basis(basis)(j)*U[m-1]*fe.grad(basis)(i,m);
              conv *= 0.5;
              break;
            default:
              std::cerr << "EqualOrderOperators::Weak::Convection: "
                        << "Unknown form " << form << std::endl;
          }
          EM((i-1)*cmp+k,(j-1)*cmp+l) += coef*conv;
        }
}


void EqualOrderOperators::Weak::Divergence(Matrix& EM, const FiniteElement& fe,
                                           double scale, int basis, int tbasis)
{
  DivGrad<1>(EM,fe,scale,basis,tbasis);
}


void EqualOrderOperators::Weak::Gradient(Matrix& EM, const FiniteElement& fe,
                                         double scale, int basis, int tbasis)
{
  DivGrad<2>(EM,fe,scale,basis,tbasis);
}


void EqualOrderOperators::Weak::Divergence(Vector& EV, const FiniteElement& fe,
                                           const Vec3& D, double scale, int basis)
{
  for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
    double div=0.0;
    for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
      div += D[k-1]*fe.grad(basis)(i,k);
    EV(i) += scale*div*fe.detJxW;
  }
}


void EqualOrderOperators::Weak::Gradient(Vector& EV, const FiniteElement& fe,
                                         double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  for (size_t i = 1; i <= fe.basis(basis).size(); ++i)
    for (size_t k = 1; k <= nsd; ++k)
      EV((i-1)*nsd+k) += scale*fe.grad(basis)(i,k)*fe.detJxW;
}


void EqualOrderOperators::Weak::Laplacian(Matrix& EM, const FiniteElement& fe,
                                          double scale, bool stress, int basis)
{
  size_t cmp = EM.rows() / fe.basis(basis).size();
  Matrix A;
  A.multiply(fe.grad(basis),fe.grad(basis),false,true);
  A *= scale*fe.detJxW;
  addComponents(EM, A, cmp, cmp, 0);
  if (stress)
    for (size_t i = 1; i <= fe.basis(basis).size(); i++)
      for (size_t j = 1; j <= fe.basis(basis).size(); j++)
        for (size_t k = 1; k <= cmp; k++)
          for (size_t l = 1; l <= cmp; l++)
            EM(cmp*(j-1)+k,cmp*(i-1)+l) += scale*fe.grad(basis)(i,k)*fe.grad(basis)(j,l)*fe.detJxW;
}


void EqualOrderOperators::Weak::LaplacianCoeff(Matrix& EM, const Matrix& K,
                                               const FiniteElement& fe,
                                               double scale, int basis)
{
  Matrix KB;
  KB.multiply(K,fe.grad(basis),false,true).multiply(scale*fe.detJxW);
  EM.multiply(fe.grad(basis),KB,false,false,true);
}


void EqualOrderOperators::Weak::Mass(Matrix& EM, const FiniteElement& fe,
                                     double scale, int basis)
{
  size_t ncmp = EM.rows()/fe.basis(basis).size();
  Matrix A;
  A.outer_product(fe.basis(basis),fe.basis(basis),false);
  A *= scale*fe.detJxW;
  addComponents(EM, A, ncmp, ncmp, 0);
}


void EqualOrderOperators::Weak::Source(Vector& EV, const FiniteElement& fe,
                                       double scale, int cmp, int basis)
{
  size_t ncmp = EV.size() / fe.basis(basis).size();
  if (cmp == 1 && ncmp == 1)
    EV.add(fe.basis(basis), scale*fe.detJxW);
  else {
    for (size_t i = 1; i <= fe.basis(basis).size(); ++i)
      for (size_t k  = (cmp == 0 ? 1: cmp);
                  k <= (cmp == 0 ? ncmp : cmp); ++k)
        EV(ncmp*(i-1)+k) += scale*fe.basis(basis)(i)*fe.detJxW;
  }
}


void EqualOrderOperators::Weak::Source(Vector& EV, const FiniteElement& fe,
                                       const Vec3& f, double scale, int basis)
{
  size_t cmp = EV.size() / fe.basis(basis).size();
  for (size_t i = 1; i <= fe.basis(basis).size(); ++i)
    for (size_t k = 1; k <= cmp; ++k)
      EV(cmp*(i-1)+k) += scale*f[k-1]*fe.basis(basis)(i)*fe.detJxW;
}


void EqualOrderOperators::Residual::Advection(Vector& EV, const FiniteElement& fe,
                                              const Vec3& AC, const Tensor& g,
                                              double scale, int basis)
{
  size_t nsd = fe.grad(basis).cols();
  for (size_t k = 1; k <= nsd; ++k) {
    double ag = g[k-1]*AC;
    for (size_t i = 1; i <= fe.basis(basis).size(); ++i)
        EV((i-1)*nsd+k) = ag*scale*fe.basis(basis)(i)*fe.detJxW;
  }
}


void EqualOrderOperators::Residual::Convection(Vector& EV, const FiniteElement& fe,
                                                const Vec3& U, const Tensor& dUdX,
                                               const Vec3& UC, double scale,
                                               WeakOperators::ConvectionForm form, int basis)
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


void EqualOrderOperators::Residual::Divergence(Vector& EV, const FiniteElement& fe,
                                                const Tensor& dUdX,
                                                double scale, size_t basis)
{
  for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
    double div=0.0;
    for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
      div += dUdX(k,k);
    EV(i) += scale*div*fe.basis(basis)(i)*fe.detJxW;
  }
}
