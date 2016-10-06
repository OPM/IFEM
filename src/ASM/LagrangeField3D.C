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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"


LagrangeField3D::LagrangeField3D (const ASMs3DLag* patch,
                                  const RealArray& v,
                                  char basis, char,
                                  const char* name) : FieldBase(name)
{
  patch->getNodalCoordinates(coord);
  patch->getSize(n1,n2,n3);
  patch->getOrder(p1,p2,p3);
  nno = n1*n2*n3;
  nelm = (n1-1)*(n2-1)*(n3-1)/(p1*p2*p3);

  // Ensure the values array has compatible length, pad with zeros if necessary
  values.resize(nno);
  RealArray::const_iterator end = v.size() > nno ? v.begin()+nno : v.end();
  std::copy(v.begin(),end,values.begin());
}


double LagrangeField3D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LagrangeField3D::valueFE (const FiniteElement& fe) const
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

  int locNode = 1;
  double value = 0.0;
  for (int k = node3; k <= node3+p3; k++)
    for (int j = node2; j <= node2+p2; j++)
      for (int i = node1; i <= node1+p1; i++, locNode++)
      {
	int node = (k-1)*n1*n2 + (j-1)*n1 + i;
	value += values(node)*N(locNode);
      }

  return value;
}


double LagrangeField3D::valueCoor (const Vec3& x) const
{
  // Not implemented yet
  return 0.0;
}


bool LagrangeField3D::gradFE (const FiniteElement& fe, Vector& grad) const
{
  grad.resize(3,true);

  Vector Vnod;
  Matrix dNdu;
  if (!Lagrange::computeBasis(Vnod,dNdu,p1,fe.xi,p2,fe.eta,p3,fe.zeta))
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
  Matrix Xnod(3,nen);

  int locNode = 1;
  for (int k = node3; k <= node3+p3; k++)
    for (int j = node2; j <= node2+p2; j++)
      for (int i = node1; i <= node1+p1; i++, locNode++)
      {
	int node = (k-1)*n1*n2 + (j-1)*n1 + i;
	Xnod.fillColumn(locNode,coord.getColumn(node));
        Vnod(locNode) = values(node);
      }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);
  return dNdX.multiply(Vnod,grad);
}


bool LagrangeField3D::gradCoor (const Vec3& x, Vector& grad) const
{
  // Not implemented yet
  return false;
}
