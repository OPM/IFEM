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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"

LagrangeFields2D::LagrangeFields2D(Matrix X, int nx, int ny, 
				   int px, int py, char* name)
  : Fields(2,name), coord(X), n1(nx), n2(ny), p1(px), p2(py) 
{
  nno = n1*n2;
  nelm = (n1-1)*(n2-1)/(p1*p2);

  // Number of fields set in fill
  nf = 0;
}


bool LagrangeFields2D::valueNode(int node, Vector& vals) const
{
  vals.resize(nf,0.0);
  for (int i = 1;i <= nf;i++)
    vals(i) = values(nf*(node-1)+i);

  return true;
}


bool LagrangeFields2D::valueFE(const FiniteElement& fe, Vector& vals) const
{
  vals.resize(nf,0.0);
  
  Vector N;
  Matrix dNdu;
  Lagrange::computeBasis(N,dNdu,p1,fe.u,p2,fe.v);

  const int nel1 = (n1-1)/p1;

  div_t divresult = div(fe.iel,nel1);
  const int iel1  = divresult.rem;
  const int iel2  = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  int node, dof;
  int locNode = 1;
  double value;
  for (int j = node2;j <= node2+p2;j++)
    for (int i = node1;i <= node1+p1;i++, locNode++) {
      node = (j-1)*n1 + i;
      value = N(locNode);
      for (int k = 1;k <= nf;k++) {
	dof = nf*(node-1) + k;
	vals(k) += values(dof)*value;
      }
    }

  return true;
}


bool LagrangeFields2D::valueCoor(const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool LagrangeFields2D::gradFE(const FiniteElement& fe, Matrix& grad) const
{
  grad.resize(nf,nsd);
  grad.fill(0.0);

  Vector N;
  Matrix dNdu;
  if (!Lagrange::computeBasis(N,dNdu,p1,fe.u,p2,fe.v))
    return false;

  const int nel1 = (n1-1)/p1;

  div_t divresult = div(fe.iel,nel1);
  const int iel1  = divresult.rem;
  const int iel2  = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;
  
  const int nen = (p1+1)*(p2+1);
  Matrix Xnod(nsd,nen);

  int node;
  int locNode = 1;
  for (int j = node2;j <= node2+p2;j++)
    for (int i = node1;i <= node1+p1;i++, locNode++) {
      node = (j-1)*n1 + i;
      Xnod.fillColumn(locNode,coord.getColumn(node));
    }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu,false);

  double value;
  int dof;
  locNode = 1;
  for (int j = node2;j <= node2+p2;j++)
    for (int i = node1;i <= node1+p1;i++, locNode++) {
      node = (j-1)*n1 + i;
      for (int k = 1;k <= nf;k++) {
	dof = nf*(node-1) + k;
	value = values(dof);
	for (int l = 1;l <= nsd;l++) 
	  grad(k,l) += value*dNdX(locNode,l);
      }
    }

  return true;
}


bool LagrangeFields2D::gradCoor(const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
