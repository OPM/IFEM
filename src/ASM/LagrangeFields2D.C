// $Id$
//==============================================================================
//!
//! \file LagrangeFields2D.C
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element vector fields in 2D.
//!
//==============================================================================

#include "LagrangeFields2D.h"
#include "ASMs2DLag.h"
#include "ItgPoint.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"


LagrangeFields2D::LagrangeFields2D (const ASMs2DLag* patch,
                                    const RealArray& v, char,
                                    const char* name)
  : Fields(name), mnpc(patch->getElmNodes(1))
{
  patch->getNodalCoordinates(coord);
  int dummy;
  patch->getOrder(p1,p2,dummy,ASM::INTEGRATION_BASIS);
  nelm = patch->getNoElms();
  nno = patch->getNoNodes();
  nf = v.size()/nno;

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nf*nno);
  RealArray::const_iterator end = v.size() > nf*nno ? v.begin()+nf*nno:v.end();
  std::copy(v.begin(),end,values.begin());
}


bool LagrangeFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+nf*(node-1));

  return true;
}


bool LagrangeFields2D::valueFE (const ItgPoint& x, Vector& vals) const
{
  vals.resize(nf,true);
  Vector N;
  std::vector<int> ip;
  if (!this->basisFE(x,N,ip))
    return false;

  for (size_t a = 0; a < ip.size(); a++)
    for (int i = 1; i <= nf; ++i)
      vals(i) += values(ip[a]*nf + i)*N[a];

  return true;
}


bool LagrangeFields2D::basisFE (const ItgPoint& x, Vector& N,
                                std::vector<int>& ip) const
{
  if (x.iel < 1 || static_cast<size_t>(x.iel) > nelm)
  {
    std::cerr <<" *** LagrangeFields2D::basisFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  if (!Lagrange::computeBasis(N,p1,x.xi,p2,x.eta))
    return false;

  ip = mnpc[x.iel-1];
  return true;
}


bool LagrangeFields2D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  grad.resize(nf,2,true);
  if (x.iel < 1 || (size_t)x.iel > nelm)
  {
    std::cerr <<" *** LagrangeFields2D::gradFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,x.xi,p2,x.eta))
    return false;

  Matrix Xnod(2,mnpc[x.iel-1].size()), Vnod(nf,mnpc[x.iel-1].size());

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
