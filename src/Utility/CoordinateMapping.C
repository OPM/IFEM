// $Id$
//==============================================================================
//!
//! \file CoordinateMapping.C
//!
//! \date May 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for coordinate mapping transformations.
//!
//==============================================================================

#include "CoordinateMapping.h"
#include "Vec3.h"
#include "Profiler.h"

#ifndef epsZ
//! \brief Zero tolerance for the Jacobian determinant.
#define epsZ 1.0e-16
#endif


Real utl::Jacobian (matrix<Real>& J, matrix<Real>& dNdX,
                    const matrix<Real>& X, const matrix<Real>& dNdu,
                    bool computeGradient)
{
  // Compute the Jacobian matrix, J = [dXdu]
  J.multiply(X,dNdu); // J = X * dNdu

  Real detJ;
  if (J.cols() == 1 && J.rows() > 1)
  {
    // Special treatment for one-parametric elements in multi-dimension space
    dNdX = dNdu;
    // Compute the length of the tangent vector {X},u
    detJ = J.getColumn(1).norm2(); // |J| = sqrt(X,u*X,u + Y,u*Y,u + Z,u*Z,u)
  }
  else if (J.cols() == 2 && J.rows() > 2)
  {
    // Special for two-parametric elements in 3-dimensional space (shell)
    dNdX = dNdu;
    // Compute the length of the normal vector {X},u x {X},v
    Vec3 a1 = J.getColumn(1);
    Vec3 a2 = J.getColumn(2);
    detJ = Vec3(a1,a2).length();
    // Store the two tangent vectors in J
    a1.normalize();
    a2.normalize();
    J.fillColumn(1,a1.ptr());
    J.fillColumn(2,a2.ptr());
  }
  else
  {
    // Compute the Jacobian determinant and inverse
    detJ = J.inverse(epsZ);

    if (computeGradient)
    {
      // Compute the first order derivatives of the basis function, w.r.t. X
      if (detJ == Real(0))
        dNdX.clear();
      else
        dNdX.multiply(dNdu,J); // dNdX = dNdu * J^-1
    }
  }

  return detJ;
}


Real utl::Jacobian (matrix<Real>& J, Vec3& t, matrix<Real>& dNdX,
                    const matrix<Real>& X, const matrix<Real>& dNdu,
                    size_t tangent)
{
  // Compute the Jacobian matrix, J = [dXdu]
  J.multiply(X,dNdu); // J = X * dNdu

  // Extract the tangent vector
  t = J.getColumn(tangent);

  // Compute the Jacobian determinant and inverse
  Real detJ = J.inverse(epsZ);

  // Compute the first order derivatives of the basis function, w.r.t. X
  if (detJ == Real(0))
    dNdX.clear();
  else
    dNdX.multiply(dNdu,J); // dNdX = dNdu * J^-1

  // Return the curve dilation (dS) in the tangent direction, vt
  return t.normalize();
}


Real utl::Jacobian (matrix<Real>& J, Vec3& n, matrix<Real>& dNdX,
                    const matrix<Real>& X, const matrix<Real>& dNdu,
                    size_t t1, size_t t2)
{
  // Compute the Jacobian matrix, J = [dXdu]
  J.multiply(X,dNdu); // J = X * dNdu

  Real dS;
  if (J.cols() == 2)
  {
    // Compute the face normal
    Vec3 v1(J.getColumn(t1));
    Vec3 v2(J.getColumn(t2));
    Vec3 v3(v1,v2); // v3 = v1 x v2
    // Compute the curve dilation (dS) and the in-plane edge normal (n)
    dS = v2.normalize(); // dA = |v2|
    v3.normalize();
    n.cross(v2,v3); // n = v2 x v3 / (|v2|*|v3|)
    if (J.rows() > 2)
    {
      // Special treatment for two-parametric elements in 3D space (shells)
      dNdX = dNdu;
      v1.normalize();
      J.fillColumn(t1,v1.ptr());
      J.fillColumn(t2,v2.ptr());
      return dS;
    }
  }
  else
  {
    // Compute the face normal (n) and surface dilation (dS)
    n.cross(J.getColumn(t1),J.getColumn(t2)); // n = t1 x t2
    dS = n.normalize(); // dS = |n|
  }

  // Compute the Jacobian inverse
  if (J.inverse(epsZ) == Real(0))
  {
    dS = Real(0);
    dNdX.clear();
  }
  else
    // Compute the first order derivatives of the basis function, w.r.t. X
    dNdX.multiply(dNdu,J); // dNdX = dNdu * J^-1

  return dS;
}


bool utl::Hessian (matrix3d<Real>& H, matrix3d<Real>& d2NdX2,
                   const matrix<Real>& Ji, const matrix<Real>& X,
                   const matrix3d<Real>& d2Ndu2, const matrix<Real>& dNdX,
                   bool geoMapping)
{
  PROFILE4("utl::Hessian");

  // Compute the Hessian matrix, H = [d2Xdu2]
  if (geoMapping && !H.multiply(X,d2Ndu2)) // H = X * d2Ndu2
    return false;
  else if (dNdX.empty())
  {
    // Probably a singular point, silently ignore
    d2NdX2.resize(0,0,0,true);
    return true;
  }
  else if (Ji.cols() <= 2 && Ji.rows() > Ji.cols())
  {
    // Special treatment for one-parametric elements in multi-dimension space
    d2NdX2 = d2Ndu2;
    return true;
  }

  // Check that the matrix dimensions are compatible
  size_t nsd = X.rows();
  if (Ji.rows() != nsd || Ji.cols() != nsd)
  {
    std::cerr <<"Hessian: Invalid dimension on Jacobian inverse, Ji("
              << Ji.rows() <<","<< Ji.cols() <<"), nsd="<< nsd << std::endl;
    return false;
  }

  // Compute the second order derivatives of the basis functions, w.r.t. X
  d2NdX2.resize(dNdX.rows(),nsd,nsd,true);
  size_t i1, i2, i3, i4, i6;
  for (size_t n = 1; n <= dNdX.rows(); n++)
    for (i1 = 1; i1 <= nsd; i1++)
      for (i2 = 1; i2 <= i1; i2++)
      {
        Real& v = d2NdX2(n,i1,i2);
        for (i3 = 1; i3 <= nsd; i3++)
          for (i4 = 1; i4 <= nsd; i4++)
          {
            Real Ji31x42 = Ji(i3,i1)*Ji(i4,i2);
            v += d2Ndu2(n,i3,i4)*Ji31x42;
            for (i6 = 1; i6 <= nsd; i6++)
              v -= dNdX(n,i6)*H(i6,i3,i4)*Ji31x42;
          }

        if (i2 < i1)
          d2NdX2(n,i2,i1) = v; // symmetry
      }

  return true;
}


void utl::getGmat (const matrix<Real>& Ji, const Real* du, matrix<Real>& G)
{
  size_t nsd = Ji.cols();
  G.resize(nsd,nsd,true);

  Real domain = pow(2.0,nsd);
  for (size_t k = 1; k <= nsd; k++)
    for (size_t l = 1; l <= nsd; l++)
    {
      Real scale = domain/(du[k-1]*du[l-1]);
      for (size_t m = 1; m <= nsd; m++)
        G(k,l) += Ji(m,k)*Ji(m,l)*scale;
    }
}
