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
#include "FiniteElement.h"
#include "Lagrange.h"
#include "CoordinateMapping.h"
#include "Vec3.h"


LagrangeFields2D::LagrangeFields2D (const ASMs2DLag* patch, const RealArray& v,
				    const char* name) : Fields(name)
{
  patch->getNodalCoordinates(coord);
  patch->getSize(n1,n2);
  patch->getOrder(p1,p2);
  nno = n1*n2;
  nelm = (n1-1)*(n2-1)/(p1*p2);
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


bool LagrangeFields2D::valueFE (const FiniteElement& fe, Vector& vals) const
{
  vals.resize(nf,true);

  Vector N;
  if (!Lagrange::computeBasis(N,p1,fe.xi,p2,fe.eta))
    return false;

  const int nel1 = (n1-1)/p1;

  div_t divresult = div(fe.iel,nel1);
  int iel1 = divresult.rem;
  int iel2 = divresult.quot;
  const int node1 = p1*iel1-1;
  const int node2 = p2*iel2-1;

  int locNode = 1;
  for (int j = node2; j <= node2+p2; j++)
    for (int i = node1; i <= node1+p1; i++, locNode++)
    {
      int dof = nf*(n1*(j-1) + i-1) + 1;
      double value = N(locNode);
      for (int k = 1; k <= nf; k++, dof++)
	vals(k) += values(dof)*value;
    }

  return true;
}


bool LagrangeFields2D::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool LagrangeFields2D::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  grad.resize(nf,2,true);

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
  Matrix Xnod(2,nen), Vnod(nf,nen);

  int locNode = 1;
  for (int j = node2; j <= node2+p2; j++)
    for (int i = node1; i <= node1+p1; i++, locNode++)
    {
      int node = (j-1)*n1 + i;
      Xnod.fillColumn(locNode,coord.getColumn(node));
      Vnod.fillColumn(locNode,values.ptr()+nf*(node-1));
    }

  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);
  grad.multiply(Vnod,dNdX);

  return true;
}


bool LagrangeFields2D::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
