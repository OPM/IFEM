// $Id$
//==============================================================================
//!
//! \file ASMs2Dmxrecovery.C
//!
//! \date Dec 13 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 2D mixed spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2Dmx.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include <array>


bool ASMs2Dmx::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                   const IntegrandBase& integrand,
                                   bool continuous) const
{
  const size_t nnod = projBasis->numCoefs_u()*projBasis->numCoefs_v();

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int p11 = projBasis->order_u();
  const int p21 = projBasis->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : nullptr;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  std::array<RealArray,2> gpar;
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg);

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPtsSf> spl0;
  std::array<std::vector<Go::BasisDerivsSf>,2> spl1;
  if (continuous) {
    projBasis->computeBasisGrid(gpar[0],gpar[1],spl1[0]);
    surf->computeBasisGrid(gpar[0],gpar[1],spl1[1]);
  } else
    projBasis->computeBasisGrid(gpar[0],gpar[1],spl0);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar.data()))
  {
    std::cerr <<" *** ASMs2Dmx::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size() << std::endl;
    return false;
  }

  double dA = 1.0;
  std::array<Vector, 2> phi;
  phi[0].resize(p11*p21);
  phi[1].resize(p1*p2);
  std::array<Matrix,2> dNdu;
  Matrix Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i2 = 0; i2 < nel2; i2++)
    for (int i1 = 0; i1 < nel1; i1++, iel++)
    {
      if (MLGE[iel] < 1) continue; // zero-area element

      if (continuous)
      {
        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,1+iel))
          return false;
        else if ((dA = 0.25*this->getParametricArea(1+iel)) < 0.0)
          return false; // topology error (probably logic error)
      }

      int ip = (i2*ng1*nel1 + i1)*ng2;
      IntVec lmnpc;
      if (projBasis != m_basis[0]) {
        lmnpc.reserve(phi[0].size());
        int vidx = (spl1[0][ip].left_idx[1]-p21+1)*projBasis->numCoefs_u();
        for (int j = 0; j < p21; ++j, vidx += projBasis->numCoefs_u())
          for (int i = 0; i < p11; ++i)
            if (continuous)
              lmnpc.push_back(spl1[0][ip].left_idx[0]-p11+1+i+vidx);
            else
              lmnpc.push_back(spl0[ip].left_idx[0]-p11+1+i+vidx);
      }
      const IntVec& mnpc = projBasis == m_basis[0] ? MNPC[iel] : lmnpc;

      // --- Integration loop over all Gauss points in each direction ----------

      Matrix eA(p11*p21, p11*p21);
      Vectors eB(sField.rows(), Vector(p11*p21));
      for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
        for (int i = 0; i < ng1; i++, ip++)
        {
          if (continuous) {
            SplineUtils::extractBasis(spl1[0][ip],phi[0],dNdu[0]);
            SplineUtils::extractBasis(spl1[1][ip],phi[1],dNdu[1]);
          }
          else
            phi[0] = spl0[ip].basisValues;

          // Compute the Jacobian inverse and derivatives
          double dJw = 1.0;
          if (continuous)
          {
            dJw = dA*wg[i]*wg[j]*utl::Jacobian(J,dNdu[1],Xnod,dNdu[1],false);
            if (dJw == 0.0) continue; // skip singular points
          }

          // Integrate the mass matrix
          eA.outer_product(phi[0], phi[0], true, dJw);

          // Integrate the rhs vector B
          for (size_t r = 1; r <= sField.rows(); r++)
            eB[r-1].add(phi[0],sField(r,ip+1)*dJw);
        }

      for (int i = 0; i < p11*p21; ++i) {
        for (int j = 0; j < p11*p21; ++j)
          A(mnpc[i]+1, mnpc[j]+1) += eA(i+1, j+1);

        int jp = mnpc[i]+1;
        for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
          B(jp) += eB[r](1+i);
      }
    }

  return true;
}
