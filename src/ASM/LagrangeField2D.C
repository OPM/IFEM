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


LagrangeField2D::LagrangeField2D (const ASMs2DLag* patch,
                                  const RealArray& v, char, char,
                                  const char* name) : FieldBase(name)
{
  patch->getNodalCoordinates(coord);
  patch->getOrder(p1,p2,n2);
  patch->getSize(n1,n2);
  nelm = (n1-1)*(n2-1)/(p1*p2);
  nno = n1*n2;

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nno);
  RealArray::const_iterator end = v.size() > nno ? v.begin()+nno : v.end();
  std::copy(v.begin(),end,values.begin());
}


double LagrangeField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LagrangeField2D::valueFE (const ItgPoint& x) const
{
  if (x.iel < 1 || (size_t)x.iel > nelm)
  {
    std::cerr <<" *** LagrangeField2D::valueFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return 0.0;
  }

  Vector N;
  Lagrange::computeBasis(N,p1,x.xi,p2,x.eta);

  div_t divresult = div(x.iel,(n1-1)/p1);
  int iel1 = divresult.rem;
  int iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  int locNode = 1;
  double value = 0.0;
  for (int j = node2; j <= node2+p2; j++)
    for (int i = node1; i <= node1+p1; i++, locNode++)
    {
      int node = (j-1)*n1 + i;
      value += values(node)*N(locNode);
    }

  return value;
}


bool LagrangeField2D::gradFE (const ItgPoint& x, Vector& grad) const
{
  grad.resize(2,true);
  if (x.iel < 1 || (size_t)x.iel > nelm)
  {
    std::cerr <<" *** LagrangeField2D::gradFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector Vnod;
  Matrix dNdu;
  if (!Lagrange::computeBasis(Vnod,dNdu,p1,x.xi,p2,x.eta))
    return false;

  div_t divresult = div(x.iel,(n1-1)/p1);
  int iel1 = divresult.rem;
  int iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  const int nen = (p1+1)*(p2+1);
  Matrix Xnod(2,nen);

  int locNode = 1;
  for (int j = node2; j <= node2+p2; j++)
    for (int i = node1; i <= node1+p1; i++, locNode++)
    {
      int node = (j-1)*n1 + i;
      Xnod.fillColumn(locNode,coord.getColumn(node));
      Vnod(locNode) = values(node);
    }

  Matrix Jac, dNdX;
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false; // Singular Jacobian

  return dNdX.multiply(Vnod,grad);
}
