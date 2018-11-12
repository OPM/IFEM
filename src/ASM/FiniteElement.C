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
