// $Id$
//==============================================================================
//!
//! \file ASMu3Dmxrecovery.C
//!
//! \date Nov 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Recovery techniques for unstructured mixed LR B-splines.
//!
//==============================================================================

#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Element.h"

#include "ASMu3Dmx.h"
#include "GlbL2projector.h"
#include "IntegrandBase.h"
#include "SparseMatrix.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#endif
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SplineUtils.h"
#include "Profiler.h"
#include <numeric>


/*!
  \brief Expands a tensor parametrization point to an unstructured one.
*/

static void expandTensorGrid (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size()*in[2].size());
  out[1].resize(in[0].size()*in[1].size()*in[2].size());
  out[2].resize(in[0].size()*in[1].size()*in[2].size());

  size_t i, j, k, ip = 0;
  for (k = 0; k < in[2].size(); k++)
    for (j = 0; j < in[1].size(); j++)
      for (i = 0; i < in[0].size(); i++, ip++) {
        out[0][ip] = in[0][i];
        out[1][ip] = in[1][j];
        out[2][ip] = in[2][k];
      }
}


bool ASMu3Dmx::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                   const IntegrandBase& integrand,
                                   bool continuous) const
{
  size_t nnod = this->getNoProjectionNodes();

  const int p1 = projBasis->order(0);
  const int p2 = projBasis->order(1);
  const int p3 = projBasis->order(2);

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const int ng3 = continuous ? nGauss : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : nullptr;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  double dV = 0.0;
  Vector phiG, phiP;
  Matrix sField, dNdu, Xnod, Jac;
  Go::BasisPts    spl0G, spl0P;
  Go::BasisDerivs spl1G, spl1P;


  // === Assembly loop over all elements in the patch ==========================

  LR::LRSplineVolume* geomBasis = m_basis[geoBasis-1].get();
  for (const LR::Element* el1 : geomBasis->getAllElements())
  {
    double uh = (el1->umin()+el1->umax())/2.0;
    double vh = (el1->vmin()+el1->vmax())/2.0;
    double wh = (el1->wmin()+el1->wmax())/2.0;
    size_t ielG = 1 + m_basis[geoBasis-1]->getElementContaining(uh,vh,wh);
    size_t ielP = 1 + projBasis->getElementContaining(uh,vh,wh);

    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,ielG))
        return false;
      else if ((dV = 0.25*this->getParametricVolume(ielG)) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    RealArray gpar[3], unstrGpar[3];
    this->getGaussPointParameters(gpar[0],0,ng1,ielG,xg);
    this->getGaussPointParameters(gpar[1],1,ng2,ielG,yg);
    this->getGaussPointParameters(gpar[2],2,ng3,ielG,zg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar, unstrGpar);

    // Evaluate the secondary solution at all integration points
    if (!this->evalSolution(sField,integrand,unstrGpar))
      return false;

    // set up basis function size (for extractBasis subroutine)
    const LR::Element* elm = projBasis->getElement(ielP-1);
    phiP.resize(elm->nBasisFunctions());
    phiG.resize(el1->nBasisFunctions());
    IntVec mnpc; mnpc.reserve(phiP.size());
    if (projBasis == m_basis[0])
      mnpc = MNPC[ielG-1];
    else for (const LR::Basisfunction* f : elm->support())
      mnpc.push_back(f->getId());

    // --- Integration loop over all Gauss points in each direction ----------

    Matrix eA(phiP.size(),phiP.size());
    Vectors eB(sField.rows(),Vector(phiP.size()));
    int ip = 0;
    for (int k = 0; k < ng3; k++)
      for (int j = 0; j < ng2; j++)
        for (int i = 0; i < ng1; i++, ip++)
        {
          if (continuous)
          {
            projBasis->computeBasis(gpar[0][i], gpar[1][j], gpar[2][k],
                                    spl1P, ielP-1);
            phiP = spl1P.basisValues;
            geomBasis->computeBasis(gpar[0][i], gpar[1][j], gpar[2][k],
                                    spl1G, ielG-1);
            SplineUtils::extractBasis(spl1G, phiG, dNdu);
          }
          else
          {
            projBasis->computeBasis(gpar[0][i], gpar[1][j], gpar[2][k],
                                    spl0P, ielP-1);
            phiP = spl0P.basisValues;
          }

          // Compute the Jacobian inverse and derivatives
          double dJw = 1.0;
          if (continuous)
          {
            dJw = dV*wg[i]*wg[j]*wg[k]*utl::Jacobian(Jac,dNdu,Xnod,dNdu,false);
            if (dJw == 0.0) continue; // skip singular points
          }

          // Integrate the mass matrix
          eA.outer_product(phiP, phiP, true, dJw);

          // Integrate the rhs vector B
          for (size_t r = 1; r <= sField.rows(); r++)
            eB[r-1].add(phiP,sField(r,ip+1)*dJw);
        }

    for (size_t i = 0; i < eA.rows(); i++)
    {
      int ip = mnpc[i]+1;
      for (size_t j = 0; j < eA.cols(); j++)
        A(ip, mnpc[j]+1) += eA(i+1,j+1);

      for (size_t k = 0; k < eB.size(); k++, ip += nnod)
        B(ip) += eB[k][i];
    }
  }

  return true;
}


bool ASMu3Dmx::globalL2projection (Matrix& sField,
                                   const IntegrandBase& integrand,
                                   bool continuous, bool /*later enforceEnds*/) const
{
  if (m_basis.empty())
    return true; // silently ignore empty patches

  PROFILE2("ASMu3Dmx::globalL2");

  // Assemble the projection matrices
  size_t nvar = std::inner_product(nb.begin(),nb.end(),nfx.begin(),0u);
  SparseMatrix* A;
  StdVector* B;
#ifdef HAS_PETSC
  if (GlbL2::MatrixType == SystemMatrix::PETSC && GlbL2::SolverParams)
  {
    A = new PETScMatrix(ProcessAdm(),*GlbL2::SolverParams,LinAlg::SYMMETRIC);
    B = new PETScVector(ProcessAdm(),nvar);
  }
  else
#endif
  {
    A = new SparseMatrix(SparseMatrix::SUPERLU);
    B = new StdVector(nvar);
  }
  A->redim(nvar,nvar);

  if (!this->assembleL2matrices(*A,*B,integrand,continuous))
  {
    delete A;
    delete B;
    return false;
  }

#if SP_DEBUG > 1
  std::cout <<"---- Matrix A -----\n"<< *A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----\n"<< *B
            <<"-------------------"<< std::endl;
#endif

  // Solve the patch-global equation system
  if (!A->solve(*B))
  {
    delete A;
    delete B;
    return false;
  }

  // Store the control-point values of the projected field
  sField.resize(*std::max_element(nfx.begin(),nfx.end()),nnod);

  size_t i, j, inod = 1, jnod = 1;
  for (size_t b = 0; b < nb.size(); b++)
    for (i = 1; i <= nb[b]; i++, inod++)
      for (j = 1; j <= nfx[b]; j++, jnod++)
        sField(j,inod) = (*B)(jnod);

#if SP_DEBUG > 1
  std::cout <<"- Solution Vector -"<< sField
            <<"-------------------"<< std::endl;
#endif
  delete A;
  delete B;
  return true;
}
