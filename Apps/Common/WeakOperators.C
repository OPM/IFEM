//==============================================================================
//!
//! \file WeakOperators.C
//!
//! \date Jul 22 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete operators.
//!
//==============================================================================

#include "WeakOperators.h"
#include "FiniteElement.h"
#include "Vec3.h"

namespace WeakOperators {
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


  //! \brief Helper applying a divergence (1) or a gradient (2) or both (0) operation
  template<int Operation>
  static void DivGrad(Matrix& EM, const FiniteElement& fe,
              double scale=1.0, size_t nf=1)
  {
    for (size_t i = 1; i <= fe.N.size();++i)
      for (size_t j = 1; j <= fe.N.size();++j)
        for (size_t k = 1; k <= fe.dNdX.cols(); ++k) {
          double div = fe.N(j)*fe.dNdX(i,k)*fe.detJxW;
          if (Operation == 2 || Operation == 0)
            EM(nf*(i-1)+k,nf*j) += -scale*div;
          if (Operation == 1 || Operation == 0)
            EM(nf*j, nf*(i-1)+k) += scale*div;
        }
  }


  void Advection(Matrix& EM, const FiniteElement& fe,
                 const Vec3& AC, double scale,
                 size_t cmp, size_t nf, size_t scmp)
  {
    Matrix C(fe.N.size(), fe.N.size());
    for (size_t i = 1; i <= fe.N.size(); ++i) {
      for (size_t j = 1; j <= fe.N.size(); ++j) {
        // Sum convection for each direction
        for (size_t k = 1; k <= fe.dNdX.cols(); ++k)
          C(i,j) += AC[k-1]*fe.dNdX(j,k);

        C(i,j) *= scale*fe.N(i)*fe.detJxW;
      }
    }
    addComponents(EM, C, cmp, nf, scmp);
  }

  void Convection(Matrix& EM, const FiniteElement& fe,
                  const Vec3& U, const Matrix& dUdX,
                  double scale, size_t cmp, size_t nf,
                  bool conservative)
  {
    if (conservative) {
      Advection(EM, fe, U, -scale, cmp, nf);
      for (size_t i = 1;i <= fe.N.size();i++)
        for (size_t j = 1;j <= fe.N.size();j++) {
          for (size_t k = 1;k <= cmp;k++) {
            for (size_t l = 1;l <= cmp;l++)
              EM((j-1)*nf+l,(i-1)*nf+k) -= scale*U[l-1]*fe.N(i)*fe.dNdX(j,k)*fe.detJxW;
          }
        }
    }
    else {
      Advection(EM, fe, U, scale, cmp, nf);
      for (size_t i = 1;i <= fe.N.size();i++)
        for (size_t j = 1;j <= fe.N.size();j++) {
          for (size_t k = 1;k <= cmp;k++) {
            for (size_t l = 1;l <= cmp;l++)
              EM((j-1)*nf+l,(i-1)*nf+k) += scale*dUdX(l,k)*fe.N(j)*fe.N(i)*fe.detJxW;
          }
        }
    }
  }


  void Divergence(Matrix& EM, const FiniteElement& fe,
                  size_t nf, double scale)
  {
    DivGrad<1>(EM,fe,scale,nf);
  }


  void PressureDiv(Matrix& EM, const FiniteElement& fe,
                   size_t nf, double scale)
  {
    DivGrad<0>(EM,fe,scale,nf);
  }


  void Gradient(Matrix& EM, const FiniteElement& fe,
                size_t nf, double scale)
  {
    DivGrad<2>(EM,fe,scale,nf);
  }


  void Divergence(Vector& EV, const FiniteElement& fe,
                  const Vec3& D, double scale, size_t cmp, size_t nf)
  {
    for (size_t i = 1; i <= fe.N.size(); ++i) {
      double div=0.0;
      for (size_t k = 1; k <= fe.dNdX.cols(); ++k)
        div += D[k-1]*fe.dNdX(i,k);
      EV((i-1)*nf+cmp) += scale*div*fe.detJxW;
    }
  }


  void Gradient(Vector& EV, const FiniteElement& fe,
                double scale, size_t nf)
  {
    for (size_t i = 1; i <= fe.N.size(); ++i)
      for (size_t k = 1; k <= fe.dNdX.cols(); ++k)
        EV((i-1)*nf+k) += scale*fe.dNdX(i,k)*fe.detJxW;
  }


  void Laplacian(Matrix& EM, const FiniteElement& fe,
                 double scale, size_t cmp, size_t nf,
                 bool stress, size_t scmp, unsigned char basis)
  {
    Matrix A;
    A.multiply(fe.grad(basis),fe.grad(basis),false,true);
    A *= scale*fe.detJxW;
    addComponents(EM, A, cmp, nf, scmp);
    if (stress)
      for (size_t i = 1; i <= fe.N.size(); i++)
        for (size_t j = 1; j <= fe.N.size(); j++)
          for (size_t k = 1; k <= cmp; k++)
            for (size_t l = 1; l <= cmp; l++)
              EM(nf*(j-1)+k+scmp,nf*(i-1)+l+scmp) += scale*fe.grad(basis)(i,k)*fe.grad(basis)(j,l)*fe.detJxW;
  }


  void LaplacianCoeff(Matrix& EM, const Matrix& K,
                      const FiniteElement& fe,
                      double scale)
  {
    Matrix KB;
    KB.multiply(K,fe.dNdX,false,true).multiply(scale*fe.detJxW);
    EM.multiply(fe.dNdX,KB,false,false,true);
  }


  void Mass(Matrix& EM, const FiniteElement& fe,
            double scale, size_t cmp, size_t nf, size_t scmp,
            unsigned char basis)
  {
    Matrix A;
    A.outer_product(fe.basis(basis),fe.basis(basis),false);
    A *= scale*fe.detJxW;
    addComponents(EM, A, cmp, nf, scmp);
  }


  void Source(Vector& EV, const FiniteElement& fe,
              double scale, size_t cmp, size_t nf, size_t scmp,
              unsigned char basis)
  {
    if (cmp == 1 && nf == 1)
      EV.add(fe.basis(basis), scale*fe.detJxW);
    else {
      for (size_t i = 1; i <= fe.N.size(); ++i)
        for (size_t k = 1; k <= cmp; ++k)
          EV(nf*(i-1)+k+scmp) += scale*fe.basis(basis)(i)*fe.detJxW;
    }
  }


  void Source(Vector& EV, const FiniteElement& fe,
              const Vec3& f, double scale, size_t cmp, size_t nf, size_t scmp)
  {
    for (size_t i=0;i<cmp;++i)
      Source(EV, fe, f[i+scmp]*scale, 1, nf, scmp+i);
  }
}
