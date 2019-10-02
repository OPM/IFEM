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
                                    const char* name) : Fields(name)
{
  patch->getNodalCoordinates(coord);
  patch->getSize(n1,n2,n3);
  patch->getOrder(p1,p2,p3);
  nno = n1*n2*n3;
  nelm = (n1-1)*(n2-1)*(n3-1)/(p1*p2*p3);
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
  if (x.iel < 1 || (size_t)x.iel > nelm)
  {
    std::cerr <<" *** LagrangeFields3D::valueFE: Element index "<< x.iel
              <<" out of range [1,"<< nelm <<"]."<< std::endl;
    return false;
  }

  Vector N;
  if (!Lagrange::computeBasis(N,p1,x.xi,p2,x.eta,p3,x.zeta))
    return false;

  const int nel1 = (n1-1)/p1;
  const int nel2 = (n2-1)/p2;

  div_t divresult = div(x.iel,nel1*nel2);
  int iel2 = divresult.rem;
  int iel3 = divresult.quot;
  divresult = div(iel2,nel1);
  int iel1 = divresult.rem;
  iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;
  const int node3 = p3*iel3-1;

  int locNode = 1;
  for (int k = node3; k <= node3+p3; k++)
    for (int j = node2; j <= node2+p2; j++)
      for (int i = node1; i <= node1+p1; i++, locNode++)
      {
	int dof = nf*(n1*(n2*(k-1) + j-1) + i-1) + 1;
	double value = N(locNode);
	for (int l = 1; l <= nf; l++, dof++)
	  vals(l) += values(dof)*value;
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

  const int nel1 = (n1-1)/p1;
  const int nel2 = (n2-1)/p2;

  div_t divresult = div(x.iel,nel1*nel2);
  int iel2 = divresult.rem;
  int iel3 = divresult.quot;
  divresult = div(iel2,nel1);
  int iel1 = divresult.rem;
  iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;
  const int node3 = p3*iel3-1;

  const int nen = (p1+1)*(p2+1)*(p3+1);
  Matrix Xnod(3,nen), Vnod(nf,nen);

  int locNode = 1;
  for (int k = node3; k <= node3+p3; k++)
    for (int j = node2; j <= node2+p2; j++)
      for (int i = node1; i <= node1+p1; i++, locNode++)
      {
	int node = (k-1)*n1*n2 + (j-1)*n1 + i;
	Xnod.fillColumn(locNode,coord.getColumn(node));
	Vnod.fillColumn(locNode,values.ptr()+nf*(node-1));
      }

  Matrix Jac, dNdX;
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false; // Singular Jacobian

  return !grad.multiply(Vnod,dNdX).empty();
}
