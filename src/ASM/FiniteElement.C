// $Id$
//==============================================================================
//!
//! \file FiniteElement.C
//!
//! \date Sep 05 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Finite element quantities at an integration point.
//!
//==============================================================================

#include "FiniteElement.h"
#include "BasisFunctionVals.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"


std::ostream& FiniteElement::write (std::ostream& os) const
{
  os <<"FiniteElement: iel="<< iel <<" iGP="<< iGP
     <<"\n               p, q, r: "<< p <<" "<< q <<" "<< r
     <<"\n               u, v, w: "<< u <<" "<< v <<" "<< w
     <<"\n               xi, eta, zeta: "<< xi <<" "<< eta <<" "<< zeta
     <<"\n               h, detJxW: "<< h <<" "<< detJxW << std::endl;
  for (size_t n = 0; n < XC.size(); n++)
    os <<"               XC_"<< n+1 <<": "<< XC[n] << std::endl;
  if (!N.empty())      os <<"N:"<< N;
  if (!dNdX.empty())   os <<"dNdX:"<< dNdX;
  if (!d2NdX2.empty()) os <<"d2NdX2: "<< d2NdX2;
  if (!d3NdX3.empty()) os <<"d3NdX3: "<< d3NdX3;
  if (!G.empty())      os <<"G:"<< G;
  if (!H.empty())      os <<"H:"<< H;
  if (!Navg.empty())   os <<"Navg:"<< Navg;
  if (!Xn.empty())     os <<"Xn:"<< Xn;
  if (!Te.isZero(0.0)) os <<"Te:\n"<< Te;
  for (size_t i = 0; i < Tn.size(); i++)
    os <<"Tn_"<< i+1 <<":\n"<< Tn[i];
  if (!P.empty())     os << "P:\n"<< P;
  if (!dPdX.empty())  os << "dPdX:\n"<< dPdX;
  return os;
}


MxFiniteElement::MxFiniteElement (const std::vector<size_t>& n, size_t ip)
  : FiniteElement(n.front(),ip)
{
  M.resize(n.size()-1);
  for (size_t b = 1; b < n.size(); b++)
    M[b-1].resize(n[b]);

  dMdX.resize(M.size());
  d2MdX2.resize(M.size());
  d3MdX3.resize(M.size());
}


std::ostream& MxFiniteElement::write (std::ostream& os) const
{
  this->FiniteElement::write(os);
  for (size_t b = 0; b < M.size(); b++)
  {
    os <<"Basis "<< b+2 <<":\n";
    if (!M[b].empty())      os <<"N:"<< M[b];
    if (!dMdX[b].empty())   os <<"dNdX:"<< dMdX[b];
    if (!d2MdX2[b].empty()) os <<"d2NdX2: "<< d2MdX2[b];
    if (!d3MdX3[b].empty()) os <<"d3NdX3: "<< d3MdX3[b];
  }
  return os;
}


/*!
  This method also calculates the first-derivatives of the basis functions
  with respect to the Cartesian coordinates, using the same geometry mapping
  for all bases.
*/

bool MxFiniteElement::Jacobian (Matrix& Jac, const Matrix& Xnod,
                                unsigned short int gBasis,
                                const BasisValuesPtrs& bfs)
{
  const bool separateGeometry = bfs.size() > this->getNoBasis();
  if (separateGeometry)
    gBasis = bfs.size();

  detJxW = utl::Jacobian(Jac, this->grad(gBasis), Xnod,
                         bfs[gBasis-1]->dNdu,
                         !separateGeometry);
  if (detJxW == 0.0) return false; // singular point

  for (size_t b = 1; b <= this->getNoBasis(); b++)
    if (b != gBasis || separateGeometry)
      this->grad(b).multiply(bfs[b-1]->dNdu, Jac);

  return true;
}


/*!
  If the input argument \a nBasis is half (or less than half) of the size of
  the internal basis function value array, it is used to flag that we are doing
  element interface terms, and the basis function values of both elements
  sharing the interface are stored internally. If \a nBasis is zero (default),
  it is reset to the size of the \a dNxdu argument, which will fit normal cases.
*/

bool MxFiniteElement::Jacobian (Matrix& Jac, Vec3& n,
                                const Matrix& Xnod,
                                unsigned short int gBasis,
                                const BasisValuesPtrs& bfs,
                                size_t t1, size_t t2, size_t nBasis,
                                const Matrix* Xnod2)
{
  if (nBasis == 0) nBasis = bfs.size();
  const bool separateGeometry = nBasis > this->getNoBasis();
  if (separateGeometry) gBasis = nBasis;

  Matrix dummy;
  if (2*nBasis <= this->getNoBasis())
  {
    // We are are doing interface terms, evaluate for the second element first
    detJxW = utl::Jacobian(Jac,n,this->grad(nBasis+gBasis),
                           Xnod2 ? *Xnod2 : Xnod,
                           bfs[nBasis+gBasis-1]->dNdu,t1,t2);

    for (size_t b = 1; b <= nBasis; ++b)
      if (b != gBasis)
        this->grad(nBasis+b).multiply(bfs[nBasis+b-1]->dNdu,Jac);
  }
  else if (separateGeometry)
    nBasis = this->getNoBasis();

  Matrix& dX = separateGeometry ? dummy : this->grad(gBasis);
  detJxW = utl::Jacobian(Jac,n,dX,Xnod,bfs[gBasis-1]->dNdu,t1,t2);

  for (size_t b = 1; b <= nBasis; ++b)
    if (b != gBasis || separateGeometry)
      this->grad(b).multiply(bfs[b-1]->dNdu,Jac);

  return detJxW != 0.0;
}


/*!
  This method also calculates the second-derivatives of the basis functions
  with respect to the Cartesian coordinates, using the same geometry mapping
  for all bases.
*/

bool MxFiniteElement::Hessian (Matrix3D& Hess, const Matrix& Jac,
                               const Matrix& Xnod,
                               unsigned short int gBasis,
                               const BasisValuesPtrs& bfs)
{
  const bool separateGeometry = bfs.size() > this->getNoBasis();
  bool ok;
  if (separateGeometry)
    ok = Hess.multiply(Xnod, bfs.back()->d2Ndu2);
  else
    ok = utl::Hessian(Hess, this->hess(gBasis), Jac, Xnod,
                      bfs[gBasis-1]->d2Ndu2,
                      this->grad(gBasis));

  for (size_t b = 1; b <= this->getNoBasis() && ok; b++)
    if (b != gBasis || separateGeometry)
      ok = utl::Hessian(Hess, this->hess(b), Jac, Xnod,
                        bfs[b-1]->d2Ndu2,
                        this->grad(b), false);

  return ok;
}


void MxFiniteElement::piolaMapping (const double detJ,
                                    const Matrix& Ji,
                                    const Matrix& Xnod,
                                    const BasisValuesPtrs& bfs)
{
  Matrix J;
  J.multiply(Xnod,bfs.back()->dNdu);
  this->piolaBasis(detJ, J);
  this->piolaGradient(detJ, J, Ji, Xnod, bfs);
}


void MxFiniteElement::piolaBasis (const double detJ, const Matrix& J)
{
  size_t np = 0, dim = J.rows();
  for (size_t b = 1; b <= dim; ++b)
    np += this->basis(b).size();

  Matrix Ntmp(dim,np);
  size_t k = 1;
  for (size_t b = 1; b <= dim; ++b)
    for (size_t i = 1; i <= this->basis(b).size(); ++i, ++k)
      Ntmp(b,k) = this->basis(b)(i);

  P.multiply(J,Ntmp,false,false,false,1.0/detJ);
}


void MxFiniteElement::piolaGradient (const double detJ,
                                     const Matrix& J, const Matrix& Ji,
                                     const Matrix& Xnod,
                                     const BasisValuesPtrs& bfs)
{
  Matrix3D H;
  H.multiply(Xnod,bfs.back()->d2Ndu2);

  size_t np = 0, dim = J.rows();
  for (size_t b = 1; b <= dim; ++b)
    np += this->basis(b).size();

  dPdX.resize(dim*dim,np);

  Matrices dJdX;
  Vector dDetdX;
  utl::JacobianGradient(Ji,H,dJdX);
  utl::detJacGradient(J,Ji,H,dDetdX);

  std::vector<Matrix> Ptmp(dim);
  for (size_t i = 0; i < dim; ++i)
  {
    Ptmp[i] = dJdX[i];
    Ptmp[i] *= 1.0 / detJ;
    Ptmp[i].add(J, -dDetdX[i] / (detJ * detJ));
  }

  size_t k = 0;
  for (size_t b = 1; b <= dim; k += this->basis(b).size(), ++b)
    for (size_t i = 1; i <= this->basis(b).size(); ++i)
    {
      Vector bf(dim);
      Matrix dBf(dim,dim);
      bf(b) = this->basis(b)(i);
      for (size_t j = 1; j <= dim; ++j)
        dBf(b,j) = bfs[b-1]->dNdu(i,j);
      for (size_t d = 1; d <= dim; ++d)
      {
        Vector tmp(dim), res(dim);
        dBf.multiply(Ji.getColumn(d), tmp);
        J.multiply(tmp, res, 1.0 / detJ);
        Ptmp[d-1].multiply(bf, res, false, 1);
        for (size_t j = 1; j <= dim; ++j)
          dPdX((d-1) * dim + j,k + i) = res(j);
      }
    }
}
