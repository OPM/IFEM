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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"


LagrangeField2D::LagrangeField2D(Matrix X, int nx, int ny, 
				 int px, int py, char* name)
  : Field(2,name), coord(X), n1(nx), n2(ny), p1(px), p2(py) 
{
  nno = n1*n2;
  nelm = (n1-1)*(n2-1)/(p1*p2);
}


double LagrangeField2D::valueNode(int node) const
{
  return values(node);
}


double LagrangeField2D::valueFE(const FiniteElement& fe) const
{
  Vector N;
  Lagrange::computeBasis(N,p1,fe.xi,p2,fe.eta);

  const int nel1 = (n1-1)/p1;

  div_t divresult = div(fe.iel,nel1);
  int iel1 = divresult.rem;
  int iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  int node, i, j;
  int locNode = 1;
  double value = 0.0;
  for (j = node2; j <= node2+p2; j++)
    for (i = node1; i <= node1+p1; i++, locNode++)
    {
      node = (j-1)*n1 + i;
      value += values(node)*N(locNode);
    }

  return value;
}


double LagrangeField2D::valueCoor(const Vec3& x) const
{
  // Not implemented yet
  return 0.0;
}


bool LagrangeField2D::gradFE(const FiniteElement& fe, Vector& grad) const
{
  grad.resize(nsd,0.0);

  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,fe.xi,p2,fe.eta))
    return false;

  const int nel1 = (n1-1)/p1;

  div_t divresult = div(fe.iel,nel1);
  int iel1 = divresult.rem;
  int iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  const int nen = (p1+1)*(p2+1);
  Matrix Xnod(nsd,nen);

  int node, i, j;
  int locNode = 1;
  for (j = node2; j <= node2+p2; j++)
    for (i = node1; i <= node1+p1; i++, locNode++)
    {
      node = (j-1)*n1 + i;
      Xnod.fillColumn(locNode,coord.getColumn(node));
    }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu,false);

  locNode = 1;
  for (j = node2; j <= node2+p2; j++)
    for (i = node1; i <= node1+p1; i++, locNode++)
    {
      node = (j-1)*n1 + i;
      grad.add(dNdX.getColumn(locNode),values(node));
    }

  return true;
}


bool LagrangeField2D::gradCoor(const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
