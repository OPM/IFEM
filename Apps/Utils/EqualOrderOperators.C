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

#include <iostream>


namespace
{
  //! \brief Helper for adding an element matrix to one component.
  //! \param[out] EM The element matrix to add to
  //! \param[in] A The scalar element matrix to add
  //! \param[in] cmp1 First components to add matrix to
  //! \param[in] cmp2 Second components to add matrix to
  //! \param[in] nf Number of components in total matrix
  void addComponent (Matrix& EM, const Matrix& A,
                     size_t cmp1, size_t cmp2, size_t nf)
  {
    for (size_t i = 1; i <= A.rows(); ++i)
      for (size_t j = 1; j <= A.cols(); ++j)
        EM(nf*(i-1)+cmp1, nf*(j-1)+cmp2) += A(i, j);
  }

  //! \brief Helper applying a divergence (1) or a gradient (2) operation.
  template<int Operation>
  void DivGrad (Matrix& EM, const FiniteElement& fe,
                double scale, int basis, int tbasis)
  {
    size_t nsd = fe.grad(basis).cols();
    for (size_t i = 1; i <= fe.basis(tbasis).size(); ++i)
      for (size_t j = 1; j <= fe.basis(basis).size(); ++j)
        for (size_t k = 1; k <= nsd; ++k) {
          double div = fe.basis(basis)(j) * fe.grad(tbasis)(i,k) * fe.detJxW;
          if (Operation == 2)
            EM((i-1)*nsd+k, j) -= scale*div;
          if (Operation == 1)
            EM(j, (i-1)*nsd+k) += scale*div;
        }
  }
}


std::pair<double,double> WeakOperators::FormConvectionWeights (ConvectionForm form)
{
  switch (form) {
    case WeakOperators::CONSERVATIVE:  return { 0.0, -1.0};
    case WeakOperators::SKEWSYMMETRIC: return { 0.5, -0.5};
    default:                           return { 1.0,  0.0};
  }
}


void EqualOrderOperators::Weak::Advection (Matrix& EM, const FiniteElement& fe,
                                           const Vec3& AC,
                                           double scale,
                                           WeakOperators::ConvectionForm form,
                                           int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();
  const Vector& N = fe.basis(basis);
  const Matrix& D = fe.grad(basis);

  // The term is made of the piece which advects the trial function and the
  // piece which advects the test function
  const auto [cw, tw] = WeakOperators::FormConvectionWeights(form);

  // One component at a time and nothing across them, so a vector problem
  // builds the one matrix they share and spreads it over the diagonal
  Matrix tmp;
  Matrix& C = cmp > 1 ? tmp : EM;
  bool addTo = cmp == 1;

  for (size_t k = 1; k <= D.cols(); ++k) {
    const double a = scale*AC[k-1]*fe.detJxW;
    if (cw != 0.0) {
      C.outer_product(N, D.getColumn(k), addTo, cw*a);
      addTo = true;
    }
    if (tw != 0.0) {
      C.outer_product(D.getColumn(k), N, addTo, tw*a);
      addTo = true;
    }
  }

  if (cmp > 1)
    for (size_t k = 1; k <= cmp; ++k)
      addComponent(EM, tmp, k, k, cmp);
}


void EqualOrderOperators::Weak::Convection (Matrix& EM, const FiniteElement& fe,
                                            const Vec3& U, const Tensor& dUdX,
                                            double scale,
                                            WeakOperators::ConvectionForm form,
                                            int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();
  const double coef = scale*fe.detJxW;

  const Vector& N = fe.basis(basis);
  const Matrix& D = fe.grad(basis);

  const auto [cw, tw] = WeakOperators::FormConvectionWeights(form);

  Matrix NN;
  if (cw != 0.0)
    NN.outer_product(N, N);

  Vector DU(N.size());
  for (size_t m = 1; m <= cmp; ++m)
    DU.add(D.getColumn(m), U[m-1]);

  Matrix C(N.size(), N.size());
  for (size_t k = 1; k <= cmp; ++k)
    for (size_t l = 1; l <= cmp; ++l) {
      C.fill(0.0);
      if (cw != 0.0) C.add(NN, cw*dUdX(k,l));
      if (tw != 0.0) C.outer_product(D.getColumn(l), N, true, tw*U[k-1]);
      if (k == l) {
        if (cw != 0.0) C.outer_product(N, DU, true, cw);
        if (tw != 0.0) C.outer_product(DU, N, true, tw);
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
  const size_t nsd = fe.grad(basis).cols();
  fe.grad(basis).multiply(Vector(D.ptr(),nsd),EV,scale*fe.detJxW,1.0);
}


void EqualOrderOperators::Weak::Laplacian (Matrix& EM, const FiniteElement& fe,
                                           double scale, bool stress, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();
  if (cmp > 1)
  {
    Matrix A;
    A.multiply(fe.grad(basis),fe.grad(basis),false,true,false,scale*fe.detJxW);
    for (size_t k = 1; k <= cmp; ++k)
      addComponent(EM, A, k, k, cmp);
  }
  else
    EM.multiply(fe.grad(basis),fe.grad(basis),false,true,true,scale*fe.detJxW);

  if (stress)
    Stress(EM, fe, scale, basis);
}


void EqualOrderOperators::Weak::Stress (Matrix& EM,
                                        const FiniteElement& fe,
                                        double scale, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();

  Matrix A;
  for (size_t k = 1; k <= cmp; k++)
    for (size_t l = 1; l <= cmp; l++) {
      A.outer_product(fe.grad(basis).getColumn(l),
                      fe.grad(basis).getColumn(k),
                      false, scale*fe.detJxW);
      addComponent(EM, A, k, l, cmp);
    }
}


void EqualOrderOperators::Weak::GradDiv (Matrix& EM,
                                         const FiniteElement& fe,
                                         double scale, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();

  Matrix A;
  for (size_t k = 1; k <= cmp; k++)
    for (size_t l = 1; l <= cmp; l++) {
      A.outer_product(fe.grad(basis).getColumn(k),
                      fe.grad(basis).getColumn(l),
                      false, scale*fe.detJxW);
      addComponent(EM, A, k, l, cmp);
    }
}


void EqualOrderOperators::Weak::LaplacianCoeff (Matrix& EM, const Matrix& K,
                                                const FiniteElement& fe,
                                                double scale, int basis)
{
  Matrix KB;
  KB.multiply(K,fe.grad(basis),false,true).multiply(scale*fe.detJxW);

  const size_t cmp = EM.rows() / fe.basis(basis).size();
  if (cmp > 1)
  {
    Matrix A;
    A.multiply(fe.grad(basis),KB,false,false);
    for (size_t k = 1; k <= cmp; ++k)
      addComponent(EM, A, k, k, cmp);
  }
  else
    EM.multiply(fe.grad(basis),KB,false,false,true);
}


void EqualOrderOperators::Weak::LaplacianStreamline (Matrix& EM,
                                                     const FiniteElement& fe,
                                                     const Vec3& AC,
                                                     double scale, int basis)
{
  const Matrix& dNdX = fe.grad(basis);
  const size_t nbf = fe.basis(basis).size();
  const size_t cmp = EM.rows() / nbf;
  const size_t nsd = dNdX.cols();

  // The streamline derivative of each test function
  Vector w(nbf);
  for (size_t k = 1; k <= nsd; ++k)
    w.add(dNdX.getColumn(k), AC[k-1]);

  for (size_t j = 1; j <= nbf; ++j)
  {
    double lap = 0.0;
    for (size_t k = 1; k <= nsd; ++k)
      lap += fe.hess(basis)(j,k,k);

    for (size_t i = 1; i <= nbf; ++i)
    {
      const double v = scale*w(i)*lap*fe.detJxW;
      for (size_t k = 1; k <= cmp; ++k)
        EM(cmp*(i-1)+k, cmp*(j-1)+k) += v;
    }
  }
}


void EqualOrderOperators::Weak::SourceStreamline (Vector& EV,
                                                  const FiniteElement& fe,
                                                  const Vec3& AC, const Vec3& f,
                                                  double scale, int basis)
{
  const Matrix& dNdX = fe.grad(basis);
  const size_t cmp = EV.size() / fe.basis(basis).size();

  Vector w(dNdX.rows());
  for (size_t k = 1; k <= dNdX.cols(); ++k)
    w.add(dNdX.getColumn(k), AC[k-1]);

  for (size_t k = 1; k <= cmp; ++k)
    EV.add(w, scale*f[k-1]*fe.detJxW, 0, 1, k-1, cmp);
}


void EqualOrderOperators::Weak::LaplacianCoeff (Matrix& EM, const Vec3& K,
                                                const FiniteElement& fe,
                                                double scale, int basis)
{
  Matrix BK(fe.grad(basis));
  BK.scale((K*scale*fe.detJxW).vec(fe.grad(basis).cols()));
  EM.multiply(BK,fe.grad(basis),false,true,true);
}


bool EqualOrderOperators::Weak::MassCoeff (Matrix& EM, const Matrix& C,
                                           const FiniteElement& fe,
                                           double scale, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();
  const double coef = scale * fe.detJxW;

  if (C.rows() < cmp || C.cols() < cmp) {
    std::cerr <<" *** EqualOrderOperators::Weak::MassCoeff: The coefficient is "
              << C.rows() <<"x"<< C.cols() <<", which does not carry the "
              << cmp <<" components of the field."<< std::endl;
    return false;
  }

  Matrix A;
  for (size_t k = 1; k <= cmp; ++k)
    for (size_t l = 1; l <= cmp; ++l)
      if (C(k,l) != 0.0) {
        A.outer_product(fe.basis(basis), fe.basis(basis), false, coef*C(k,l));
        if (cmp > 1)
          addComponent(EM, A, k, l, cmp);
        else
          EM.add(A);
      }

  return true;
}


void EqualOrderOperators::Weak::Mass (Matrix& EM, const FiniteElement& fe,
                                      double scale, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();
  if (cmp > 1)
  {
    Matrix A;
    A.outer_product(fe.basis(basis),fe.basis(basis),false,scale*fe.detJxW);
    for (size_t k = 1; k <= cmp; ++k)
      addComponent(EM, A, k, k, cmp);
  }
  else
    EM.outer_product(fe.basis(basis),fe.basis(basis),true,scale*fe.detJxW);
}


void EqualOrderOperators::Weak::ItgConstraint (Matrix& EM, const FiniteElement& fe,
                                               double scale, int basis)
{
  const size_t cmp = EM.rows() / fe.basis(basis).size();

  Matrix A(fe.basis(basis).size(), 1);
  A = fe.basis(basis);
  A *= scale*fe.detJxW;
  for (size_t k = 1; k <= cmp; ++k)
    addComponent(EM, A, k, k, cmp);
}


void EqualOrderOperators::Weak::Source (Vector& EV,
                                        const FiniteElement& fe,
                                        double scale, int cmp, int basis)
{
  const size_t ncmp = EV.size() / fe.basis(basis).size();

  if (cmp == 1 && ncmp == 1)
    EV.add(fe.basis(basis), scale*fe.detJxW);
  else if (cmp > 0)
    EV.add(fe.basis(basis), scale*fe.detJxW, 0, 1, cmp-1, ncmp);
  else for (size_t k = 0; k < ncmp; k++)
    EV.add(fe.basis(basis), scale*fe.detJxW, 0, 1, k, ncmp);
}


void EqualOrderOperators::Weak::Source (Vector& EV, const FiniteElement& fe,
                                        const Vec3& f, double scale, int basis)
{
  const size_t cmp = EV.size() / fe.basis(basis).size();

  for (size_t k = 1; k <= cmp; ++k)
    EV.add(fe.basis(basis), scale*f[k-1]*fe.detJxW, 0, 1, k-1, cmp);
}


void EqualOrderOperators::Residual::Advection (Vector& EV,
                                               const FiniteElement& fe,
                                               const Vec3& AC, const Tensor& g,
                                               double scale, int basis)
{
  const size_t nsd = fe.grad(basis).cols();

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
  const size_t cmp = EV.size() / fe.basis(basis).size();
  const double coef = scale * fe.detJxW;

  const Vector& N = fe.basis(basis);
  const Matrix& D = fe.grad(basis);

  // The term is made of the piece which advects the velocity and the piece
  // which advects the test function
  const auto [cw, tw] = WeakOperators::FormConvectionWeights(form);

  for (size_t i = 1; i <= N.size(); i++)
    for (size_t k = 1; k <= cmp; k++)
      for (size_t l = 1; l <= cmp; l++) {
        double conv = 0.0;
        if (cw != 0.0) conv -= cw*UC[l-1]*dUdX(k,l)*N(i);
        if (tw != 0.0) conv -= tw*U[k-1]*UC[l-1]*D(i,l);
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
  const size_t nsd = fe.grad(basis).cols();

  for (size_t k = 1; k <= nsd; ++k)
    EV.add(fe.grad(basis).getColumn(k), scale*fe.detJxW, 0, 1, k-1, nsd);
}


void EqualOrderOperators::Residual::Laplacian (Vector& EV,
                                               const FiniteElement& fe,
                                               const Vec3& dUdX,
                                               double scale, int basis)
{
  const size_t nsd = fe.grad(basis).cols();

  fe.grad(basis).multiply(Vector(dUdX.ptr(),nsd),EV,-scale*fe.detJxW,1.0);
}


void EqualOrderOperators::Residual::Laplacian (Vector& EV,
                                               const FiniteElement& fe,
                                               const Tensor& dUdX,
                                               double scale,
                                               bool stress, int basis)
{
  const size_t nsd = fe.grad(1).cols();

  Tensor dUdXT(dUdX,true);
  for (size_t k = 1; k <= nsd; ++k) {
    Vector diff;
    fe.grad(basis).multiply(Vector(dUdXT[k-1].ptr(), nsd), diff);
    if (stress)
      fe.grad(basis).multiply(Vector(dUdX[k-1].ptr(), nsd), diff, false, 1);
    EV.add(diff, -scale*fe.detJxW, 0, 1, k-1, nsd);
  }
}
