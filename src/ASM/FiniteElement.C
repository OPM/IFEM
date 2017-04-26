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


std::ostream& FiniteElement::write (std::ostream& os) const
{
  os <<"FiniteElement: iel="<< iel <<" iGP="<< iGP <<" p="<< p
     <<"\n               u, v, w: "<< u <<" "<< v <<" "<< w
     <<"\n               xi, eta, zeta: "<< xi <<" "<< eta <<" "<< zeta
     <<"\n               h, detJxW: "<< h <<" "<< detJxW << std::endl;
  if (!N.empty())      os <<"N:"<< N;
  if (!dNdX.empty())   os <<"dNdX:"<< dNdX;
  if (!d2NdX2.empty()) os <<"d2NdX2:"<< d2NdX2;
  if (!G.empty())      os <<"G:"<< G;
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
  Nx.resize(n.size()-1);
  for (size_t b = 1; b < n.size(); b++)
    Nx[b-1].resize(n[b]);

  dNxdX.resize(Nx.size());
  d2NxdX2.resize(Nx.size());
}


std::ostream& MxFiniteElement::write (std::ostream& os) const
{
  this->FiniteElement::write(os);
  for (size_t b = 0; b < Nx.size(); b++)
  {
    os <<"Basis "<< b+2 <<":\n";
    if (!Nx[b].empty())      os <<"N:"<< Nx[b];
    if (!dNxdX[b].empty())   os <<"dNdX:"<< dNxdX[b];
    if (!d2NxdX2[b].empty()) os <<"d2NdX2:"<< d2NxdX2[b];
  }
  return os;
}
