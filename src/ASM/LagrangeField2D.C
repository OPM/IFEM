// $Id$
//==============================================================================
//!
//! \file LagrangeField2D.C
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element scalar field in 2D.
//!
//==============================================================================

#include "LagrangeField2D.h"
#include "ASMs2DLag.h"
#include "ItgPoint.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"

#include <numeric>


LagrangeField2D::LagrangeField2D (const ASMs2DLag* patch,
                                  const RealArray& v, char, char cmp,
                                  const char* name)
  : FieldBase(name), mnpc(patch->getElmNodes(1))
{
  patch->getNodalCoordinates(coord);
  int nf;
  patch->getOrder(p1,p2,nf);

  nelm = patch->getNoElms();
  nno = patch->getNoNodes();

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nno);

  nf = patch->getNoFields(1);
  size_t ndof = nf > 1 && cmp > 0 ? nf*nno : nno;
  RealArray::const_iterator vit = v.begin();
  RealArray::const_iterator end = v.size() > ndof ? vit+ndof : v.end();
  if (nf == 1 || cmp == 0)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nf)
      values[i] = *(vit+cmp-1);
}


double LagrangeField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LagrangeField2D::valueFE (const ItgPoint& x) const
{
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeField2D::valueFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return 0.0;
  }

  Vector N;
  Lagrange::computeBasis(N,p1,x.xi,p2,x.eta);

  int locNode = 0;
  return std::accumulate(mnpc[x.iel-1].begin(), mnpc[x.iel-1].end(), 0.0,
                         [&locNode, &v = values, &N]
                         (const double acc, const int node)
                         { return acc + v[node]*N(++locNode); });
}


bool LagrangeField2D::gradFE (const ItgPoint& x, Vector& grad) const
{
  grad.resize(2,true);
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeField2D::gradFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector Vnod;
  Matrix dNdu;
  if (!Lagrange::computeBasis(Vnod,dNdu,p1,x.xi,p2,x.eta))
    return false;

  Matrix Xnod(2,mnpc[x.iel-1].size());

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
