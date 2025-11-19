// $Id$
//==============================================================================
//!
//! \file LagrangeFields3D.C
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element vector fields in 3D.
//!
//==============================================================================

#include "LagrangeFields3D.h"
#include "ASMs3DLag.h"
#include "ItgPoint.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"


LagrangeFields3D::LagrangeFields3D (const ASMs3DLag* patch,
                                    const RealArray& v, char,
                                    const char* name)
  : Fields(name), mnpc(patch->getElmNodes(1))
{
  patch->getNodalCoordinates(coord);
  patch->getOrder(p1,p2,p3);
  nno = patch->getNoNodes();
  nelm = patch->getNoElms();
  nf = v.size()/nno;

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nf*nno);
  RealArray::const_iterator end = v.size() > nf*nno ? v.begin()+nf*nno:v.end();
  std::copy(v.begin(),end,values.begin());
}


bool LagrangeFields3D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+nf*(node-1));

  return true;
}


bool LagrangeFields3D::valueFE (const ItgPoint& x, Vector& vals) const
{
  vals.resize(nf,true);
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeFields3D::valueFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector N;
  if (!Lagrange::computeBasis(N,p1,x.xi,p2,x.eta,p3,x.zeta))
    return false;

  int locNode = 1;
  for (const int node : mnpc[x.iel-1]) {
    for (int i = 1; i <= nf; ++i)
      vals(i) +=  values(node*nf + i)* N(locNode);
    ++locNode;
  }

  return true;
}


bool LagrangeFields3D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  grad.resize(nf,3,true);
  if (x.iel < 1 || (size_t)x.iel > nelm)
  {
    std::cerr <<" *** LagrangeFields3D::gradFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,x.xi,p2,x.eta,p3,x.zeta))
    return false;

  Matrix Xnod(3,mnpc[x.iel-1].size()), Vnod(nf,mnpc[x.iel-1].size());

  int locNode = 1;
  for (const int node : mnpc[x.iel-1]) {
    for (int i = 1; i <= nf; ++i) {
      Xnod.fillColumn(locNode,coord.getColumn(node+1));
      Vnod.fillColumn(locNode,values.ptr()+nf*node);
    }
    ++locNode;
  }

  Matrix Jac, dNdX;
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false; // Singular Jacobian

  return !grad.multiply(Vnod,dNdX).empty(); // grad = Vnod * dNdX
}
