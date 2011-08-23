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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"


LagrangeField3D::LagrangeField3D(Matrix X, int nx, int ny, int nz, 
				 int px, int py, int pz, char* name)
  : Field(2,name), coord(X), n1(nx), n2(ny), n3(nz),
    p1(px), p2(py), p3(pz) 
{
  nno = n1*n2*n3;
  nelm = (n1-1)*(n2-1)*(n3-1)/(p1*p2*p3);
}


double LagrangeField3D::valueNode(int node) const
{
  return values(node);
}


double LagrangeField3D::valueFE(const FiniteElement& fe) const
{
  Vector N;
  Lagrange::computeBasis(N,p1,fe.xi,p2,fe.eta,p3,fe.zeta);

  const int nel1 = (n1-1)/p1;
  const int nel2 = (n2-1)/p2;

  div_t divresult = div(fe.iel,nel1*nel2);
  int iel2 = divresult.rem;
  int iel3 = divresult.quot;
  divresult = div(iel2,nel1);
  int iel1 = divresult.rem;
  iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;
  const int node3 = p3*iel3-1;

  int node, i, j, k;
  int locNode = 1;
  double value = 0.0;
  for (k = node3; k <= node3+p3; k++)
    for (j = node2; j <= node2+p2; j++)
      for (i = node1; i <= node1+p1; i++, locNode++)
      {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	value += values(node)*N(locNode);
      }

  return value;
}


double LagrangeField3D::valueCoor(const Vec3& x) const
{
  // Not implemented yet
  return 0.0;
}


bool LagrangeField3D::gradFE(const FiniteElement& fe, Vector& grad) const
{
  grad.resize(nsd,0.0);

  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta))
    return false;

  const int nel1 = (n1-1)/p1;
  const int nel2 = (n2-1)/p2;

  div_t divresult = div(fe.iel,nel1*nel2);
  int iel2 = divresult.rem;
  int iel3 = divresult.quot;
  divresult = div(iel2,nel1);
  int iel1 = divresult.rem;
  iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;
  const int node3 = p3*iel3-1;

  const int nen = (p1+1)*(p2+1)*(p3+1);
  Matrix Xnod(nsd,nen);

  int node, i, j, k;
  int locNode = 1;
  for (k = node3; k <= node3+p3; k++)
    for (j = node2; j <= node2+p2; j++)
      for (i = node1; i <= node1+p1; i++, locNode++) {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	Xnod.fillColumn(locNode,coord.getColumn(node));
      }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu,false);

  locNode = 1;
  for (k = node3; k <= node3+p3; k++)
    for (j = node2; j <= node2+p2; j++)
      for (i = node1;i <= node1+p1;i++, locNode++) {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	grad.add(dNdX.getColumn(locNode),values(node));
      }

  return true;
}


bool LagrangeField3D::gradCoor(const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
