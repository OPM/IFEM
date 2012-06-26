// $Id$
//==============================================================================
//!
//! \file MortarContact.C
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for Mortar-based contact analysis.
//!
//==============================================================================

#include "MortarContact.h"
#include "RigidBody.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "SAM.h"
#include "Vec3Oper.h"


MortarMats::MortarMats (const SAM& _sam, int nMast, usint n) : sam(_sam), nsd(n)
{
  int nnod = sam.getNoNodes();
  phiA.resize(nsd*nnod,nnod-nMast);
  gNA.resize(nnod);
  AA.resize(nnod);

  contact = false;
}


void MortarMats::initialize (bool)
{
  phiA.init();
  gNA.fill(0.0);
  AA.fill(0.0);
}


bool MortarMats::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (!elMat || elMat->A.size() != 1 || elMat->b.size() != 2)
    return false;

#if INT_DEBUG > 2
  std::cout <<"MortarMats: AA for element "<< elmId << elMat->b[0]
	    <<"MortarMats: gNA for element "<< elmId << elMat->b[1]
	    <<"MortarMats: phiA for element "<< elmId << elMat->A[0];
#endif

  std::vector<int> mnpc;
  if (!sam.getElmNodes(mnpc,elmId))
    return false;

  size_t i, j;
  for (j = 0; j < elMat->b[0].size() && j < mnpc.size(); j++)
    AA(mnpc[j]) += elMat->b[0][j];
  for (j = 0; j < elMat->b[1].size() && j < mnpc.size(); j++)
    gNA(mnpc[j]) += elMat->b[1][j];

  const Matrix& eM = elMat->A.front();
  for (j = 0; j < eM.cols() && j < mnpc.size(); j++)
    for (i = 0; nsd*(i+1) <= eM.rows() && i < mnpc.size(); i++)
      for (usint d = 1; d <= nsd; d++)
	phiA(nsd*(mnpc[i]-1)+d,mnpc[j]) += eM(nsd*i+d,j+1);

  return true;
}


bool MortarMats::finalize (bool)
{
  contact = false;
  for (size_t j = 1; j <= AA.size(); j++)
    if (AA(j) > 0.0)
    {
      for (size_t i = 1; i <= phiA.rows(); i++)
	if (phiA(i,j) != 0.0)
	  phiA(i,j) /= AA(j);

      gNA(j) /= AA(j);
      if (gNA(j) > 0.0)
	AA(j) = 0.0; // no contact yet
      else
	contact = true;
    }

#if INT_DEBUG > 2
  std::cout <<"MortarMats: Global AA:"<< AA
	    <<"MortarMats: Global gNA:"<< gNA
	    <<"MortarMats: Global phiA:"<< phiA;
#endif
  return true;
}


MortarContact::MortarContact (RigidBody* mst, GlobalIntegral* gi, usint nsd)
  : myInt(gi), master(mst)
{
  npv = nsd;         // Number of primary unknowns per node
  primsol.resize(1); // Only the current solution is needed
}


LocalIntegral* MortarContact::getLocalIntegral (size_t nen, size_t, bool) const
{
  size_t nedof = npv*(nen+master->getNoNodes());
  ElmMats* elm = new ElmMats;
  elm->withLHS = true;
  elm->resize(1,2);
  elm->A.front().resize(nedof,nen);
  elm->b.front().resize(nen);
  elm->b.back().resize(nen);
  return elm;
}


bool MortarContact::evalBou (LocalIntegral& elmInt,
			     const FiniteElement& fe,
			     const Vec3& X, const Vec3& sNorm) const
{
  // Evaluate the Gauss point coordinates in current configuration
  Vec3 Xd, mNorm;
  for (usint d = 0; d < npv; d++)
    Xd[d] = X[d] + elmInt.vec.front().dot(fe.N,d,npv);

  // Evaluate the gap function at current point
  Vector Nm;
  double gap = master->evalGap(Xd,mNorm,Nm);
  if (mNorm.isZero())
    mNorm = -fe.detJxW*sNorm;
  else
    mNorm *= fe.detJxW;

#if INT_DEBUG > 3
  std::cout <<"MortarContact::gap("<< Xd <<") = "<< gap << std::endl;
#endif

  // Accumulate the weighted gap matrices
  ElmMats& eMat = static_cast<ElmMats&>(elmInt);
  Matrix& phiA  = eMat.A.front();
  size_t nSlave = fe.N.size();
  for (size_t j = 1; j < nSlave; j++)
    for (size_t i = 1; i <= nSlave + Nm.size(); i++)
      for (usint k = 1; k <= npv; k++)
        if (i <= nSlave)
          phiA((i-1)*npv+k,j) += mNorm[k-1]*fe.N(i)*fe.N(j);
        else
          phiA((i-1)*npv+k,j) -= mNorm[k-1]*Nm(i-nSlave)*fe.N(j);

  eMat.b[0].add(fe.N,fe.detJxW);
  eMat.b[1].add(fe.N,gap*fe.detJxW);
  return true;
}


MortarPenalty::MortarPenalty (RigidBody* mst, const MortarMats& mats, usint nsd)
  : MortarContact(mst,NULL,nsd), mortar(mats)
{
}


LocalIntegral* MortarPenalty::getLocalIntegral (size_t nen, size_t, bool) const
{
  size_t nedof = npv*(nen+master->getNoNodes());
  ElmMats* elm = new ElmMats;
  elm->withLHS = true;
  elm->resize(1,1);
  elm->A.front().resize(nedof,nedof);
  elm->b.front().resize(nedof);
  return elm;
}


bool MortarPenalty::initElementBou (const std::vector<int>& MNPC,
				    LocalIntegral& elmInt)
{
  if (!this->IntegrandBase::initElementBou(MNPC,elmInt))
    return false;

  size_t nMast = master->getNoNodes();
  size_t nen   = MNPC.size() - nMast;
  size_t nedof = npv*(nen+nMast);
  phiA.resize(nedof,nen);
  gNA.resize(nen);

  for (size_t j = 0; j < nen; j++)
  {
    gNA(j+1) = mortar.weightedGap(MNPC[j]+1);
    for (size_t i = 0; npv*(i+1) <= nedof; i++)
      for (usint d = 1; d <= npv; d++)
	phiA(npv*i+d,j+1) = mortar.phi(npv*MNPC[i]+d,MNPC[j]+1);
  }

#if INT_DEBUG > 2
  std::cout <<"MortarPenalty::gNA"<< gNA
	    <<"MortarPenalty::phiA"<< phiA;
#endif
  return true;
}


bool MortarPenalty::evalBou (LocalIntegral& elmInt,
			     const FiniteElement& fe,
			     const Vec3& X, const Vec3& sNorm) const
{
  size_t i, j, idof, jdof;
  usint  k, l;

  // Evaluate the Gauss point coordinates in current configuration
  Vec3 Xd, mNorm(sNorm);
  for (usint k = 0; k < npv; k++)
    Xd[k] = X[k] + elmInt.vec.front().dot(fe.N,k,npv);

  // Evaluate the geometric stiffness of the contact spring
  Matrix Tg;
  Vector Nm;
  master->geometricStiffness(Xd,mNorm,Tg,Nm);

  // Evaluate p_Nint and phi_Nint at current integration point
  double epsJxW = master->eps * fe.detJxW;
  size_t nSlave = fe.N.size();
  size_t nelnod = nSlave + Nm.size();
  double p_Nint = fe.N.dot(gNA) * epsJxW;
  Matrix phi_Nint(npv,nelnod);
  for (j = jdof = 1; j <= nelnod; j++)
    for (k = 1; k <= npv; k++, jdof++)
      phi_Nint(k,j) = fe.N.dot(phiA.getRow(jdof)) * epsJxW;

#if INT_DEBUG > 3
  std::cout <<"MortarPenalty::p_Nint = "<< p_Nint
	    <<"\nMortarPenalty::phi_Nint:"<< phi_Nint;
#endif

  // Accumulate the displacement residuals

  Vector& R = static_cast<ElmMats&>(elmInt).b.front();
  for (j = jdof = 1; j <= nelnod; j++)
    for (l = 0; l < npv; l++, jdof++)
      R(jdof) += p_Nint*mNorm[l] * (j <= nSlave ? -fe.N(j) : Nm(j-nSlave));

  // Accumulate the tangent stiffness matrix

  Matrix& Kt = static_cast<ElmMats&>(elmInt).A.front();
  for (j = jdof = 1; j <= nelnod; j++)
    for (l = 1; l <= npv; l++, jdof++)
      for (i = idof = 1; i <= nelnod; i++)
      {
	double phiw = (i <= nSlave ? fe.N(i) : -Nm(i-nSlave)) * phi_Nint(l,j);
	for (k = 0; k < npv; k++, idof++)
	  Kt(idof,jdof) += phiw*mNorm[k];
      }

  if (Tg.normInf() > 1.0e-12)
  {
    // Add geometric stiffness terms
    Tg *= p_Nint;
    for (j = jdof = 1; j <= nelnod; j++)
    {
      double Nj = j <= nSlave ? fe.N(j) : Nm(j-nSlave);
      for (l = 1; l <= npv; l++, jdof++)
	for (i = idof = 1; i <= nelnod; i++)
	{
	  double Ni = i <= nSlave ? fe.N(i) : Nm(i-nSlave);
	  for (k = 1; k <= npv; k++, idof++)
	    Kt(idof,jdof) += Ni*Nj * Tg(k,l);
	}
    }
  }

#if INT_DEBUG > 3
  std::cout <<"MortarPenalty::R"<< R
	    <<"MortarPenalty::Kt"<< Kt;
#endif
  return true;
}
