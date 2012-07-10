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
#include "SAM.h"
#include "Vec3Oper.h"
#include "Utilities.h"


MortarMats::MortarMats (const SAM& _sam, int nMast, usint n) : sam(_sam), nsd(n)
{
  int nnod = sam.getNoNodes();
  phiA.resize(nsd*nnod,nnod-nMast);
  gNA.resize(nnod);
  AA.resize(nnod);
  nNoEl.resize(nnod);
  contact = false;
}


void MortarMats::initialize (bool)
{
  phiA.init();
  gNA.fill(0.0);
  AA.fill(0.0);
  std::fill(nNoEl.begin(),nNoEl.end(),0);
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
    if (mnpc[j] > 0)
    {
      AA(mnpc[j]) += elMat->b[0][j];
      gNA(mnpc[j]) += elMat->b[1][j];
      nNoEl[mnpc[j]-1] ++;
    }

  const Matrix& eM = elMat->A.front();
  for (j = 0; j < eM.cols() && j < mnpc.size(); j++)
    if (mnpc[j] > 0)
      for (i = 0; nsd*(i+1) <= eM.rows() && i < mnpc.size(); i++)
        if (mnpc[i] > 0)
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
      if (gNA(j) <= 0.0)
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


bool MortarContact::initElementBou (const std::vector<int>& MNPC,
				    LocalIntegral& elmInt)
{
  // Extract the current displacements for this element
  size_t nen = static_cast<ElmMats&>(elmInt).b.front().size();
  std::vector<int> MNPCu(MNPC.begin(),MNPC.begin()+nen);
  return this->IntegrandBase::initElementBou(MNPCu,elmInt);
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


void MortarContact::accResAndTangent (Vector& R, Matrix& Kt,
				      Vec3& mNorm, Vector& Nm,
				      const Vector& Ns, const Vec3& X,
				      const Vector& lambda, const Matrix& phi,
				      double detJxW, double epsJxW) const
{
  if (epsJxW == 0.0) epsJxW = detJxW;

  size_t i, j, idof, jdof;
  usint  k, l;

  // Evaluate the geometric stiffness of the contact spring (Tg)
  // the normal vector of the master surface (mNorm)
  // and the master node interpolation functions (Nm)
  Matrix Tg;
  master->geometricStiffness(X,mNorm,Tg,Nm);

  // Evaluate p_Nint and phi_Nint at current integration point
  size_t nSlave = Ns.size();
  size_t nelnod = nSlave + Nm.size();
  double p_Nint = Ns.dot(lambda) * detJxW;
  Matrix phi_Nint(npv,nelnod);
  for (j = jdof = 1; j <= nelnod; j++)
    for (k = 1; k <= npv; k++, jdof++)
      phi_Nint(k,j) = Ns.dot(phi.getRow(jdof)) * epsJxW;

#if INT_DEBUG > 3
  std::cout <<"MortarContact::p_Nint = "<< p_Nint
	    <<"\nMortarContact::phi_Nint:"<< phi_Nint;
#endif

  // Accumulate the residuals

  for (j = jdof = 1; j <= nelnod; j++)
    for (l = 0; l < npv; l++, jdof++)
      R(jdof) += p_Nint*mNorm[l] * (j <= nSlave ? -Ns(j) : Nm(j-nSlave));

  // Accumulate the tangent stiffness matrix

  for (j = jdof = 1; j <= nelnod; j++)
    for (l = 1; l <= npv; l++, jdof++)
      for (i = idof = 1; i <= nelnod; i++)
      {
	double phiw = (i <= nSlave ? Ns(i) : -Nm(i-nSlave)) * phi_Nint(l,j);
	for (k = 0; k < npv; k++, idof++)
	  Kt(idof,jdof) += phiw*mNorm[k];
      }

  if (Tg.normInf() > 1.0e-12)
  {
    // Add geometric stiffness terms
    Tg *= p_Nint;
    for (j = jdof = 1; j <= nelnod; j++)
    {
      double Nj = j <= nSlave ? Ns(j) : Nm(j-nSlave);
      for (l = 1; l <= npv; l++, jdof++)
	for (i = idof = 1; i <= nelnod; i++)
	{
	  double Ni = i <= nSlave ? Ns(i) : Nm(i-nSlave);
	  for (k = 1; k <= npv; k++, idof++)
	    Kt(idof,jdof) += Ni*Nj * Tg(k,l);
	}
    }
  }

#if INT_DEBUG > 3
  std::cout <<"MortarContact::R"<< R
	    <<"MortarContact::Kt"<< Kt;
#endif
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
  elm->resize(2,2);
  elm->A.front().resize(nedof,nedof);
  elm->b.front().resize(nedof);
  elm->A.back().resize(nedof,nen);
  elm->b.back().resize(nen);
  return elm;
}


bool MortarPenalty::initElementBou (const std::vector<int>& MNPC,
				    LocalIntegral& elmInt)
{
  if (!this->IntegrandBase::initElementBou(MNPC,elmInt))
    return false;

  Vector& gNA  = static_cast<ElmMats&>(elmInt).b.back();
  Matrix& phiA = static_cast<ElmMats&>(elmInt).A.back();

  for (size_t j = 1; j <= gNA.size() && j <= phiA.cols(); j++)
    if (MNPC[j-1] >= 0)
    {
      gNA(j) = mortar.weightedGap(MNPC[j-1]+1);
      if (gNA(j) > 0.0) // we have separation at this point
	gNA(j) = 0.0;
      else for (size_t i = 0; npv*(i+1) <= phiA.rows(); i++)
        if (MNPC[i] >= 0)
          for (usint d = 1; d <= npv; d++)
            phiA(npv*i+d,j) = mortar.phi(npv*MNPC[i]+d,MNPC[j-1]+1);
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
  // Evaluate the Gauss point coordinates in current configuration
  Vec3 Xd, mNorm(sNorm);
  for (usint k = 0; k < npv; k++)
    Xd[k] = X[k] + elmInt.vec.front().dot(fe.N,k,npv);

  // Accumulate displacement residual and associated tangent stiffness

  Vector Nm;
  ElmMats& elm = static_cast<ElmMats&>(elmInt);
  this->accResAndTangent(elm.b.front(),elm.A.front(),mNorm,Nm,fe.N,Xd,
			 elm.b.back(),elm.A.back(),master->eps*fe.detJxW);

  return true;
}


MortarAugmentedLag::MortarAugmentedLag (RigidBody* mst, const MortarMats& mats,
					usint nsd)
  : MortarContact(mst,NULL,nsd), mortar(mats)
{
}


LocalIntegral* MortarAugmentedLag::getLocalIntegral (size_t nen,
                                                     size_t, bool) const
{
  return new ALElmMats(npv*(nen+master->getNoNodes()));
}


bool MortarAugmentedLag::initElementBou (const std::vector<int>& MNPC,
					 LocalIntegral& elmInt)
{
  ALElmMats& elm = static_cast<ALElmMats&>(elmInt);
  size_t nedof = elm.b.front().size();
  size_t nen = nedof/npv - master->getNoNodes();

  // Extract the current displacements for this element
  std::vector<int> MNPCu(MNPC.begin(),MNPC.begin()+nen);
  if (!this->IntegrandBase::initElementBou(MNPCu,elmInt))
    return false;

  // Extract the current Lagrange multipliers for this element,. Assume one for
  // each boundary node (which should be the only nodes with positive index).
  size_t i, j, ip;
  RealArray lambda;
  for (i = ip = 0; i < nen; i++)
    if (MNPC[i] >= 0)
    {
      int ilag = MNPC[nedof/npv+(ip++)]; // Global Lagrange multiplier index
      lambda.push_back(primsol.front()[ilag]);
    }

  Matrix& Kul  = elm.A[1];
  Matrix& Kll  = elm.A[2];
  Matrix& phiA = elm.A[3];
  Vector& Rl   = elm.b[1];
  Vector& LNA  = elm.b[2];
  size_t  nlag = lambda.size();

  Kul.resize(nedof,nlag);
  Kll.resize(nlag,nlag);
  Rl.resize(nlag);

  LNA.resize(nen);
  phiA.resize(nedof,nen);
  elm.iLag.resize(nen,0);

  int lnod = 0;
  for (ip = 0, j = 1; j <= nen; j++)
    if ((lnod = MNPC[j-1]+1) > 0)
    {
      double AA = mortar.weightedArea(lnod);
      double gNA = mortar.weightedGap(lnod);
      double mult = 1.0 / mortar.getALconnectivity(lnod);

      // Check contact status at this point
      double lamA = lambda[ip++];
      if ((LNA(j) = lamA + master->eps * gNA) <= 0.0)
	elm.iLag[j-1] = ip; // we have contact, store index for Kul matrix

      // AL multiplier residual
      Rl(j) = (elm.iLag[j-1] ? gNA : -lamA/master->eps) * AA *mult;

      if (elm.iLag[j-1])
	// AL multiplier tangent stiffness
	Kll(j,j) = master->eps * AA * mult;
      else // we have separation at this point
	LNA(j) = 0.0;

      for (i = 0; npv*(i+1) <= nedof; i++)
        if (MNPC[i] >= 0)
          for (usint d = 1; d <= npv; d++)
            phiA(npv*i+d,j) = mortar.phi(npv*MNPC[i]+d,lnod);
    }

#if INT_DEBUG > 2
  std::cout <<"MortarAugmentedLag::LNA"<< LNA
	    <<"MortarAugmentedLag::phiA"<< phiA
	    <<"MortarAugmentedLag::Rl"<< Rl
	    <<"MortarAugmentedLag::Kll"<< Kll;
#endif
  return true;
}


bool MortarAugmentedLag::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& sNorm) const
{
  // Evaluate the Gauss point coordinates in current configuration
  usint k;
  Vec3 Xd, mNorm(sNorm);
  for (k = 0; k < npv; k++)
    Xd[k] = X[k] + elmInt.vec.front().dot(fe.N,k,npv);

  // Accumulate displacement residual and associated tangent stiffness matrix

  Vector Nm;
  ALElmMats& elm = static_cast<ALElmMats&>(elmInt);
  this->accResAndTangent(elm.b.front(),elm.A.front(),mNorm,Nm,fe.N,Xd,
			 elm.b[2],elm.A[3],fe.detJxW,master->eps*fe.detJxW);

  // Accumulate the tangent stiffness matrix for Augmented Lagrange multipliers

  size_t nSlave = fe.N.size();
  size_t nMast = Nm.size();
  size_t i, j, idof, jdof;
  Matrix& Kul = elm.A[1];
  for (j = 1; j <= nSlave; j++)
    if ((jdof = elm.iLag[j-1]))
    {
      // Slave surface nodes
      for (i = idof = 1; i <= nSlave; i++)
	for (k = 0; k < npv; k++, idof++)
	  Kul(idof,jdof) += mNorm[k] * fe.N(i) * fe.N(j) * fe.detJxW;

      // Master surface nodes
      for (i = 1; i <= nMast; i++)
	for (k = 0; k < npv; k++, idof++)
	  Kul(idof,jdof) -= mNorm[k] * Nm(i) * fe.N(j) * fe.detJxW;
    }

#if INT_DEBUG > 3
  std::cout <<"MortarAugmentedLag::Kul"<< Kul;
#endif
  return true;
}


MortarAugmentedLag::ALElmMats::ALElmMats (size_t nedof)
{
  withLHS = true;
  this->resize(5,4);
  A.front().resize(nedof,nedof);
  b.front().resize(nedof);
}


const Matrix& MortarAugmentedLag::ALElmMats::getNewtonMatrix () const
{
  Matrix& K = const_cast<Matrix&>(A.back());

  size_t ndof = A[0].cols();
  size_t nlag = A[1].cols();
  K.resize(ndof+nlag,ndof+nlag);

  size_t i, j;
  for (i = 1; i <= ndof; i++)
  {
    // Insert the Kuu matrix
    for (j = 1; j <= ndof; j++)
      K(i,j) = A[0](i,j);

    // Insert the Kul matrix symmetrically
    for (j = 1; j <= nlag; j++)
    {
      K(i,ndof+j) = A[1](i,j);
      K(ndof+j,i) = A[1](i,j);
    }
  }

  // Insert the Kll matrix
  for (i = 1; i <= nlag; i++)
    for (j = 1; j <= nlag; j++)
      K(ndof+i,ndof+j) = A[2](i,j);

#if INT_DEBUG > 0
  std::cout <<"\nMortarAugmentedLag::ALElmMats:Newton matrix ="<< A.back();
#endif
  return A.back();
}


const Vector& MortarAugmentedLag::ALElmMats::getRHSVector () const
{
  Vector& R = const_cast<Vector&>(b.back());

  R.clear();
  R.insert(R.end(),b[0].begin(),b[0].end()); // Ru
  R.insert(R.end(),b[1].begin(),b[1].end()); // Rl

#if INT_DEBUG > 0
  std::cout <<"\nMortarAugmentedLag::ALElmMats:RHS vector ="<< b.back();
#endif
  return b.back();
}
