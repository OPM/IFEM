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
    Vec3 dXdu(J.getColumn(1));
    detJ = dXdu.length(); // |J| = sqrt(X,u*X,u + Y,u*Y,u + Z,u*Z,u)
  }
  else if (J.cols() == 2 && J.rows() > 2)
  {
    // Special for two-parametric elements in 3-dimensional space (shell)
    dNdX = dNdu;
    // Compute the length of the normal vector {X},u x {X},v
    detJ = Vec3(J.getColumn(1),J.getColumn(2)).length();
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
    d2NdX2.clear();
    return true;
  }

  size_t nsd = X.rows();
  if (Ji.cols() < nsd)
  {
    // Special treatment for one-parametric elements in multi-dimension space
    // as well as two-parametric elements in 3D space (shells)
    d2NdX2 = d2Ndu2;
    return true;
  }

  // Check that the matrix dimensions are compatible
  if (Ji.cols() != nsd || Ji.rows() != nsd)
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


void utl::Hessian (const matrix3d<Real>& Hess, matrix<Real>& H)
{
  size_t n = Hess.dim(2);
  if (n == Hess.dim(3))
  {
    H.resize(Hess.dim(1),n*(n+1)/2);
    for (size_t i = 1; i <= n; i++)
      H.fillColumn(i,Hess.getColumn(i,i));
    for (size_t i = 2; i <= n; i++)
      H.fillColumn(n+i-1,Hess.getColumn(1,i));
    if (n > 2)
      H.fillColumn(n+n,Hess.getColumn(2,3));
  }
  else
    H.clear();
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


bool utl::Hessian2 (matrix4d<Real>& d3NdX3,
                    const matrix<Real>& Ji, const matrix4d<Real>& d3Ndu3)
{
  PROFILE4("utl::Hessian2");

  size_t nsd = Ji.rows();
  if (Ji.cols() < nsd)
  {
    // Special treatment for one-parametric elements in multi-dimension space
    // as well as two-parametric elements in 3D space (shells)
    d3NdX3 = d3Ndu3;
    return true;
  }

  // Check that the matrix dimensions are compatible
  if (Ji.cols() != nsd)
  {
    std::cerr <<"Hessian2: Invalid dimension on Jacobian inverse, Ji("
              << Ji.rows() <<","<< Ji.cols() <<")"<< std::endl;
    return false;
  }

  // Compute the third order derivatives of the basis functions, w.r.t. X
  d3NdX3.resize(d3Ndu3.dim(1),nsd,nsd,nsd,true);
  if (nsd == 1)
  {
    for (size_t i = 1; i <= d3Ndu3.dim(1); i++)
      d3NdX3(i,1,1,1) = d3Ndu3(i,1,1,1)*Ji(1,1)*Ji(1,1)*Ji(1,1);
    return true;
  }
  else if (nsd != 2)
  {
    std::cerr <<"Hessian2: Not implemented for nsd="<< nsd << std::endl;
    return false;
  }

  for (size_t i = 1; i <= d3Ndu3.dim(1); i++)
  {
    d3NdX3(i,1,1,1) =    d3Ndu3(i,1,1,1)*Ji(1,1)*Ji(1,1)*Ji(1,1) +
                         d3Ndu3(i,1,2,2)*Ji(2,1)*Ji(2,1)*Ji(1,1) +
                       2*d3Ndu3(i,1,1,2)*Ji(1,1)*Ji(1,1)*Ji(2,1) +
                         d3Ndu3(i,2,2,2)*Ji(2,1)*Ji(2,1)*Ji(2,1) +
                         d3Ndu3(i,1,1,2)*Ji(1,1)*Ji(1,1)*Ji(2,1) +
                       2*d3Ndu3(i,1,2,2)*Ji(1,1)*Ji(2,1)*Ji(2,1);

    d3NdX3(i,2,2,2) =    d3Ndu3(i,1,1,1)*Ji(1,2)*Ji(1,2)*Ji(1,2) +
                         d3Ndu3(i,1,2,2)*Ji(2,2)*Ji(2,2)*Ji(1,2) +
                       2*d3Ndu3(i,1,1,2)*Ji(1,2)*Ji(1,2)*Ji(2,2) +
                         d3Ndu3(i,1,1,2)*Ji(1,2)*Ji(1,2)*Ji(2,2) +
                         d3Ndu3(i,2,2,2)*Ji(2,2)*Ji(2,2)*Ji(2,2) +
                       2*d3Ndu3(i,1,2,2)*Ji(1,2)*Ji(2,2)*Ji(2,2);

    d3NdX3(i,1,1,2) =    d3Ndu3(i,1,1,1)*Ji(1,1)*Ji(1,1)*Ji(1,2) +
                         d3Ndu3(i,1,2,2)*Ji(2,1)*Ji(2,1)*Ji(1,2) +
                       2*d3Ndu3(i,1,1,2)*Ji(1,1)*Ji(1,2)*Ji(2,1) +
                         d3Ndu3(i,1,1,2)*Ji(1,1)*Ji(1,1)*Ji(2,2) +
                         d3Ndu3(i,2,2,2)*Ji(2,1)*Ji(2,1)*Ji(2,2) +
                       2*d3Ndu3(i,1,2,2)*Ji(1,1)*Ji(2,1)*Ji(2,2);
    d3NdX3(i,2,1,1) = d3NdX3(i,1,2,1) = d3NdX3(i,1,1,2);

    d3NdX3(i,1,2,2) =   d3Ndu3(i,1,1,1)*Ji(1,1)*Ji(1,2)*Ji(1,2) +
                        d3Ndu3(i,1,2,2)*Ji(1,1)*Ji(2,2)*Ji(2,2) +
                      2*d3Ndu3(i,1,1,2)*Ji(1,1)*Ji(1,2)*Ji(2,2) +
                        d3Ndu3(i,1,1,2)*Ji(1,2)*Ji(1,2)*Ji(2,1) +
                        d3Ndu3(i,2,2,2)*Ji(2,1)*Ji(2,2)*Ji(2,2) +
                      2*d3Ndu3(i,1,2,2)*Ji(1,2)*Ji(2,1)*Ji(2,2);
    d3NdX3(i,2,2,1) = d3NdX3(i,2,1,2) = d3NdX3(i,1,2,2);
  }

  return true;
}


void utl::JacobianGradient (const matrix<Real>& dudX,
                            const matrix3d<Real>& d2Xdu2,
                            std::vector<matrix<Real>>& dJdX)
{
  size_t i, j, k, dim = dudX.rows();

  std::vector<matrix<Real>> dJdu(dim);
  for (k = 1; k <= dim; ++k)
  {
    matrix<Real>& dJ = dJdu[k-1];
    dJ.resize(dim, dim);
    for (i = 1; i <= dim; ++i)
      for (j = 1; j <= dim; ++j)
        dJ(i,j) = d2Xdu2(i,k,j);
  }

  dJdX.resize(dim);
  for (j = 1; j <= dim; ++j)
  {
    matrix<Real>& dJ = dJdX[j-1];
    dJ.resize(dim, dim);
    for (i = 0; i < dim; ++i)
      dJ.add(dJdu[i], dudX(i+1,j));
  }
}


void utl::detJacGradient (const matrix<Real>& J, const matrix<Real>& Ji,
                          const matrix3d<Real>& H, std::vector<Real>& dDetdX)
{
  size_t dim = J.rows();

  std::vector<Real> dDetdu(dim,Real(0));

  if (dim == 2)
  {
    dDetdu[0] =   H(1,1,1)*J(2,2) + J(1,1)*H(2,2,1)
                - H(1,2,1)*J(2,1) - J(1,2)*H(2,1,1);

    dDetdu[1] =   H(1,1,2)*J(2,2) + J(1,1)*H(2,2,2)
                - H(1,2,2)*J(2,1) - J(1,2)*H(2,1,2);
  }
  else if (dim == 3)
  {
    auto det = [&J](size_t i)
    {
      const size_t i1 = (i == 1 ? 2 : 1);
      const size_t i2 = (i == 3 ? 2 : 3);
      return J(2,i1) * J(3,i2) - J(2,i2) * J(3,i1);
    };

    auto detD = [&J,&H](size_t i, size_t d)
    {
      const size_t i1 = (i == 1 ? 2 : 1);
      const size_t i2 = (i == 3 ? 2 : 3);
      return   H(2,i1,d) * J(3,i2) + J(2,i1) * H(3,i2,d)
             - H(2,i2,d) * J(3,i1) - J(2,i2) * H(3,i1,d);
    };

    for (size_t i = 1; i <= 3; ++i)
      dDetdu[i-1] =    H(1,1,i) * det(1) + J(1,1) * detD(1, i)
                    - (H(1,2,i) * det(2) + J(1,2) * detD(2, i))
                    +  H(1,3,i) * det(3) + J(1,3) * detD(3, i);
  }

  dDetdX.resize(dim,Real(0));
  Ji.multiply(dDetdu, dDetdX, true);
}
