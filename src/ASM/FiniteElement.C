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
                                const std::vector<Matrix>& dNxdu,
                                unsigned short int gBasis)
{
  detJxW = utl::Jacobian(Jac,
                         gBasis > 1 ? dMdX[gBasis-2] : dNdX,
                         Xnod, dNxdu[gBasis-1]);
  if (detJxW == 0.0) return false; // singular point

  if (gBasis > 1)
    dNdX.multiply(dNxdu.front(),Jac);

  for (size_t basis = 2; basis <= dNxdu.size(); basis++)
    if (basis != gBasis)
      dMdX[basis-2].multiply(dNxdu[basis-1],Jac);

  return true;
}


/*!
  This method also calculates the second-derivatives of the basis functions
  with respect to the Cartesian coordinates, using the same geometry mapping
  for all bases.
*/

bool MxFiniteElement::Hessian (Matrix3D& Hess, const Matrix& Jac,
                               const Matrix& Xnod,
                               const std::vector<Matrix3D>& d2Nxdu2,
                               unsigned short int gBasis)
{
  bool ok = utl::Hessian(Hess,
                         gBasis > 1 ? d2MdX2[gBasis-2] : d2NdX2,
                         Jac, Xnod, d2Nxdu2[gBasis-1],
                         gBasis > 1 ? dMdX[gBasis-2] : dNdX);

  if (ok && gBasis > 1)
    ok = utl::Hessian(Hess, d2NdX2, Jac, Xnod,
                      d2Nxdu2.front(), dNdX, false);

  for (size_t basis = 2; basis <= d2Nxdu2.size() && ok; basis++)
    if (basis != gBasis)
      ok = utl::Hessian(Hess, d2MdX2[basis-2], Jac, Xnod,
                        d2Nxdu2[basis-1], dMdX[basis-2], false);

  return ok;
}
