// $Id$
//==============================================================================
//!
//! \file LagrangeField3D.C
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element scalar field in 3D.
//!
//==============================================================================

#include "LagrangeField3D.h"
#include "ASMs3DLag.h"
#include "ItgPoint.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include <numeric>


LagrangeField3D::LagrangeField3D (const ASMs3DLag* patch,
                                  const RealArray& v, char, char cmp,
                                  const char* name)
  : FieldBase(name), mnpc(patch->getElmNodes(1))
{
  patch->getNodalCoordinates(coord);
  patch->getOrder(p1,p2,p3);
  nno = patch->getNoNodes();
  nelm = patch->getNoElms();

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nno);
  const int nf = patch->getNoFields(1);
  const size_t ndof = nf > 1 && cmp > 0 ? nf*nno : nno;
  RealArray::const_iterator vit = v.begin();
  RealArray::const_iterator end = v.size() > ndof ? vit+ndof : v.end();
  if (nf == 1 || cmp == 0)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nf)
      values[i] = *(vit+cmp-1);
}


double LagrangeField3D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LagrangeField3D::valueFE (const ItgPoint& x) const
{
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeField3D::valueFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return 0.0;
  }

  Vector N;
  Lagrange::computeBasis(N,p1,x.xi,p2,x.eta,p3,x.zeta);

  int locNode = 0;
  return std::accumulate(mnpc[x.iel-1].begin(), mnpc[x.iel-1].end(), 0.0,
                         [&locNode, &v = values, &N]
                         (const double acc, const int node)
                         { return acc + v[node]*N(++locNode); });
}


bool LagrangeField3D::gradFE (const ItgPoint& x, Vector& grad) const
{
  grad.resize(3,true);
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeField3D::gradFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector Vnod;
  Matrix dNdu;
  if (!Lagrange::computeBasis(Vnod,dNdu,p1,x.xi,p2,x.eta,p3,x.zeta))
    return false;

  Matrix Xnod(3,mnpc[x.iel-1].size());

  int locNode = 0;
  for (const int node : mnpc[x.iel-1]) {
    Xnod.fillColumn(++locNode,coord.getColumn(node+1));
    Vnod(locNode) = values[node];
  }

  Matrix Jac, dNdX;
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false; // Singular Jacobian

  return dNdX.multiply(Vnod,grad,true);
}
