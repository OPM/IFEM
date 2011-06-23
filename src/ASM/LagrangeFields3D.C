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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"


LagrangeFields3D::LagrangeFields3D(Matrix X, int nx, int ny, int nz, 
				   int px, int py, int pz, char* name)
  : Fields(3,name), coord(X), n1(nx), n2(ny), n3(nz),
    p1(px), p2(py), p3(pz) 
{
  nno = n1*n2*n3;
  nelm = (n1-1)*(n2-1)*(n3-1)/(p1*p2*p3);

  // Number of fields set in fill
  nf = 0;
}


bool LagrangeFields3D::valueNode(int node, Vector& vals) const
{
  vals.resize(nf,0.0);
  for (int i = 1;i <= nf;i++)
    vals(i) = values(nf*(node-1)+i);

  return true;
}


bool LagrangeFields3D::valueFE(const FiniteElement& fe, Vector& vals) const
{
  vals.resize(nf,0.0);
  
  Vector N;
  Matrix dNdu;
  Lagrange::computeBasis(N,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta);

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

  int node, dof;
  int locNode = 1;
  double value;
  for (int k = node3;k <= node3+p3;k++)
    for (int j = node2;j <= node2+p2;j++)
      for (int i = node1;i <= node1+p1;i++, locNode++) {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	value = N(locNode);
	for (int l = 1;l <= nf;l++) {
	  dof = nf*(node-1)+l;
	  vals(l) += values(dof)*value;
	}
      }
  
  return true;
}


bool LagrangeFields3D::valueCoor(const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool LagrangeFields3D::gradFE(const FiniteElement& fe, Matrix& grad) const
{
  grad.resize(nf,nsd);
  grad.fill(0.0);

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

  int node;
  int locNode = 1;
  for (int k = node3;k <= node3+p3;k++) 
    for (int j = node2;j <= node2+p2;j++) 
      for (int i = node1;i <= node1+p1;i++, locNode++) {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	Xnod.fillColumn(locNode,coord.getColumn(node));
      }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu,false);

  int dof;
  double value;
  locNode = 1;
  for (int k = node3;k <= node3+p3;k++) 
    for (int j = node2;j <= node2+p2;j++)
      for (int i = node1;i <= node1+p1;i++, locNode++) {
	node = (k-1)*n1*n2 + (j-1)*n1 + i;
	for (int m = 1;m <= nf;m++) {
	  dof = nf*(node-1) + m;
	  value = values(dof);
	  for (int n = 1;n <= nsd;n++)
	    grad(m,n) += value*dNdX(locNode,n);
	}
      }

  return true;
}


bool LagrangeFields3D::gradCoor(const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
