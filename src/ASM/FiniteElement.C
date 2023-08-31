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
                                const std::vector<const BasisFunctionVals*>* bf,
                                const std::vector<Matrix>* dNxdu)
{
  size_t nBasis = bf ? bf->size() : dNxdu->size();
  const bool separateGeometry = nBasis > this->getNoBasis();
  if (separateGeometry)
    gBasis = nBasis;

  detJxW = utl::Jacobian(Jac, this->grad(gBasis), Xnod,
                         bf ? (*bf)[gBasis-1]->dNdu : (*dNxdu)[gBasis-1],
                         !separateGeometry);
  if (detJxW == 0.0) return false; // singular point

  for (size_t b = 1; b <= this->getNoBasis(); b++)
    if (b != gBasis || separateGeometry)
      this->grad(b).multiply(bf ? (*bf)[b-1]->dNdu : (*dNxdu)[b-1], Jac);

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
                                const std::vector<Matrix>& dNxdu,
                                size_t t1, size_t t2, size_t nBasis,
                                const Matrix* Xnod2)
{
  if (nBasis == 0) nBasis = dNxdu.size();
  const bool separateGeometry = nBasis > this->getNoBasis();
  if (separateGeometry) gBasis = nBasis;

  Matrix dummy;
  if (2*nBasis <= this->getNoBasis())
  {
    // We are are doing interface terms, evaluate for the second element first
    detJxW = utl::Jacobian(Jac,n,this->grad(nBasis+gBasis),
                           Xnod2 ? *Xnod2 : Xnod,
                           dNxdu[nBasis+gBasis-1],t1,t2);

    for (size_t b = 1; b <= nBasis; ++b)
      if (b != gBasis)
        this->grad(nBasis+b).multiply(dNxdu[nBasis+b-1],Jac);
  }
  else if (separateGeometry)
    nBasis = this->getNoBasis();

  Matrix& dX = separateGeometry ? dummy : this->grad(gBasis);
  detJxW = utl::Jacobian(Jac,n,dX,Xnod,dNxdu[gBasis-1],t1,t2);

  for (size_t b = 1; b <= nBasis; ++b)
    if (b != gBasis || separateGeometry)
      this->grad(b).multiply(dNxdu[b-1],Jac);

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
                               const std::vector<const BasisFunctionVals*>* bf,
                               const std::vector<Matrix3D>* d2Nxdu2)
{
  size_t nBasis = bf ? bf->size() : d2Nxdu2->size();
  const bool separateGeometry = nBasis > this->getNoBasis();
  bool ok;
  if (separateGeometry)
    ok = Hess.multiply(Xnod, bf ? bf->back()->d2Ndu2 : d2Nxdu2->back());
  else
    ok = utl::Hessian(Hess, this->hess(gBasis), Jac, Xnod,
                      bf ? (*bf)[gBasis-1]->d2Ndu2 : (*d2Nxdu2)[gBasis-1],
                      this->grad(gBasis));

  for (size_t b = 1; b <= this->getNoBasis() && ok; b++)
    if (b != gBasis || separateGeometry)
      ok = utl::Hessian(Hess, this->hess(b), Jac, Xnod,
                        bf ? (*bf)[b-1]->d2Ndu2 : (*d2Nxdu2)[b-1],
                        this->grad(b), false);

  return ok;
}
