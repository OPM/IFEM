// $Id$
//==============================================================================
//!
//! \file ASMs3DLagrecovery.C
//!
//! \date Jul 7 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 3D FE models.
//!
//==============================================================================

#include "ASMs3DLag.h"
#include "ItgPoint.h"
#include "Field.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "Lagrange.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include <array>


/*
bool ASMs3DLag::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                    const IntegrandBase& integrand,
                                    bool continuous) const
{
  const size_t nnod = this->getNoNodes(1);

  const int nel1 = (nx-1)/(p1-1);
  const int nel2 = (ny-1)/(p2-1);
  const int nel3 = (nz-1)/(p3-1);

  int pmax = p1 > p2 ? p1 : p2;
  if (pmax < p3) pmax = p3;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? this->getNoGaussPt(pmax,true) : p1 - 1;
  const int ng2 = continuous ? ng1 : p2 - 1;
  const int ng3 = continuous ? ng2 : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng1) : 0;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  std::array<RealArray,3> gpar;
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg);
  gpar[2] = this->getGaussPointParameters(gp,2,ng3,zg);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar.data()))
  {
    std::cerr <<" *** ASMs3DLag::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size()*gpar[2].size()
              << std::endl;
    return false;
  }

  double dV = 1.0;
  Vector phi(p1*p2*p3);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i3 = 0; i3 < nel3; i3++)
    for (int i2 = 0; i2 < nel2; i2++)
      for (int i1 = 0; i1 < nel1; i1++, iel++)
      {
        if (MLGE[iel] < 1) continue; // zero-volume element

        if (continuous)
        {
          // Set up control point (nodal) coordinates for current element
          if (!this->getElementCoordinates(Xnod,1+iel))
            return false;
          else if ((dV = 0.125*this->getParametricVolume(1+iel)) < 0.0)
            return false; // topology error (probably logic error)
        }

        int ip = ((i3*ng2*nel2 + i2)*ng1*nel1 + i1)*ng3;
        const IntVec& mnpc = MNPC[iel];

        // --- Integration loop over all Gauss points in each direction --------

        Matrix eA(p1*p2*p3, p1*p2*p3);
        Vectors eB(sField.rows(), Vector(p1*p2*p3));
        for (int k = 0; k < ng3; k++)
          for (int j = 0; j < ng2; j++)
            for (int i = 0; i < ng1; i++, ip++)
            {
              // Compute basis function derivatives at current point
              // using tensor product of one-dimensional Lagrange polynomials
              if (continuous) {
                if (!Lagrange::computeBasis(phi, dNdu,
                                            p1,xg[i],p2,xg[j],p3,xg[k]))
                  return false;
              } else {
                if (!Lagrange::computeBasis(phi,
                                            p1,xg[i],p2,xg[j],p3,xg[k]))
                  return false;
              }

              // Compute the Jacobian inverse and derivatives
              double dJw = dV;
              if (continuous)
              {
                dJw *= wg[i]*wg[j]*wg[k]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
                if (dJw == 0.0) continue; // skip singular points
              }

              // Integrate the mass matrix
              eA.outer_product(phi, phi, true, dJw);

              // Integrate the rhs vector B
              for (size_t r = 1; r <= sField.rows(); r++)
                eB[r-1].add(phi,sField(r,ip+1)*dJw);
            }

        for (int i = 0; i < p1*p2*p3; ++i) {
          for (int j = 0; j < p1*p2*p3; ++j)
            A(mnpc[i]+1, mnpc[j]+1) += eA(i+1, j+1);

          int jp = mnpc[i]+1;
          for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
            B(jp) += eB[r](1+i);
        }
      }

  return true;
}
*/
