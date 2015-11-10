// $Id$
//==============================================================================
//!
//! \file ASMu2Dmxrecovery.C
//!
//! \date May 11 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Recovery techniques for unstructured mixed LR B-splines.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"

#include "ASMu2Dmx.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "IntegrandBase.h"


static void expandTensorGrid (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size());
  out[1].resize(in[0].size()*in[1].size());

  size_t i, j, ip = 0;
  for (j = 0; j < in[1].size(); j++)
    for (i = 0; i < in[0].size(); i++, ip++) {
      out[0][ip] = in[0][i];
      out[1][ip] = in[1][j];
    }
}


bool ASMu2Dmx::globalL2projection (Matrix& sField,
				   const IntegrandBase& integrand,
				   bool continuous) const
{
  if (!m_basis[0] || !m_basis[1]) return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::globalL2");

  const int p1 = m_basis[0]->order(0);
  const int p2 = m_basis[0]->order(1);

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : NULL;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;


  // Set up the projection matrices
  const size_t nnod = std::inner_product(nb.begin(), nb.end(), nfx.begin(), 0);
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod);
  A.redim(nnod,nnod);

  double dA = 0.0;
  Vectors phi(m_basis.size());
  std::vector<Matrix> dNdu(m_basis.size());
  Matrix Xnod, Jac;
  std::vector<Go::BasisDerivsSf> spl1(m_basis.size());
  std::vector<Go::BasisPtsSf> spl0(m_basis.size());


  // === Assembly loop over all elements in the patch ==========================

  std::vector<LR::Element*>::iterator el1 = m_basis[geoBasis-1]->elementBegin();
  for (int iel = 1; el1 != m_basis[geoBasis-1]->elementEnd(); ++el1, ++iel)
  {
    double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
    double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;
    std::vector<size_t> els;
    std::vector<size_t> elem_sizes;
    for (size_t i=0; i < m_basis.size(); ++i) {
      els.push_back(m_basis[i]->getElementContaining(uh, vh)+1);
      elem_sizes.push_back((*(m_basis[i]->elementBegin()+els.back()-1))->nBasisFunctions());
    }

    int geoEl = els[geoBasis-1];

    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,geoEl))
        return false;
      else if ((dA = 0.25*this->getParametricArea(geoEl)) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    RealArray gpar[2], unstrGpar[2];
    this->getGaussPointParameters(gpar[0],0,ng1,geoEl,xg);
    this->getGaussPointParameters(gpar[1],1,ng2,geoEl,yg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar, unstrGpar);

    // Evaluate the secondary solution at all integration points
    if (!this->evalSolution(sField,integrand,unstrGpar))
      return false;

    // set up basis function size (for extractBasis subroutine)
    for (size_t b = 0; b < m_basis.size(); ++b)
      phi[b].resize((*(m_basis[b]->elementBegin()+els[b]-1))->nBasisFunctions());

    // --- Integration loop over all Gauss points in each direction ----------
    int ip = 0;
    for (int j = 0; j < ng2; j++)
      for (int i = 0; i < ng1; i++, ip++)
      {
        if (continuous)
        {
          for (size_t b = 0; b < m_basis.size(); ++b) {
            m_basis[b]->computeBasis(gpar[0][i],gpar[1][j],spl0[b],els[b]-1);
            SplineUtils::extractBasis(spl1[b],phi[b],dNdu[b]);
          }
        }
        else
        {
          for (size_t b = 0; b < m_basis.size(); ++b) {
            m_basis[b]->computeBasis(gpar[0][i],gpar[1][j],spl0[b],els[b]-1);
            phi[b] = spl0[b].basisValues;
          }
        }

        // Compute the Jacobian inverse and derivatives
        double dJw = 1.0;
        if (continuous)
        {
          dJw = dA*wg[i]*wg[j]*utl::Jacobian(Jac,dNdu[geoBasis-1],Xnod,dNdu[geoBasis-1],false);
          if (dJw == 0.0) continue; // skip singular points
        }

        // Integrate the linear system A*x=B
        size_t el_ofs = 0;
        size_t eq_ofs = 0;
        size_t nod_ofs = 1;
        for (size_t b = 0; b < m_basis.size(); ++b) {
          for (size_t ii = 0; ii < phi[b].size(); ii++)
          {
            int inod = MNPC[iel-1][ii+el_ofs]+1;
            for (size_t jj = 0; jj < phi[b].size(); jj++)
            {
              int jnod = MNPC[iel-1][jj+el_ofs]+1;
              for (size_t k=1;k<=nfx[b];++k) {
                A((inod-nod_ofs)*nfx[b]+k+eq_ofs,(jnod-nod_ofs)*nfx[b]+k+eq_ofs) += phi[b][ii]*phi[b][jj]*dJw;
	   }
            }
            for (size_t k=1;k<=nfx[b];++k)
              B((inod-nod_ofs)*nfx[b]+k+eq_ofs) += phi[b][ii]*sField(k,ip+1)*dJw;
          }
          el_ofs += elem_sizes[b];
          eq_ofs += nb[b]*nfx[b];
	  nod_ofs += nb[b];
        }
      }
  }

#if SP_DEBUG > 2
  std::cout <<"---- Matrix A -----"<< A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----"<< B
            <<"-------------------"<< std::endl;
#endif

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the control-point values of the projected field
  sField.resize(*std::max_element(nfx.begin(), nfx.end()),
                std::accumulate(nb.begin(), nb.end(), 0));

  size_t eq_ofs = 0;
  size_t ofs = 0;
  for (size_t b = 0; b < m_basis.size(); ++b) {
    for (size_t i = 1; i <= nb[b]; i++)
      for (size_t j = 1; j <= nfx[b]; j++)
        sField(j,i+ofs) = B((i-1)*nfx[b]+j+eq_ofs);
    eq_ofs += nb[b]*nfx[b];
    ofs += nb[b];
  }

  return true;
}
