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
#include "TimeDomain.h"
#include "SAM.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include <cfloat>


MortarMats::MortarMats (const SAM& _sam, int nMast, usint n) : sam(_sam), nsd(n)
{
  int nNod = sam.getNoNodes();
  int nSlv = nNod - nMast;
  phiA.resize(nsd*nNod,nSlv);
  gNA.resize(nSlv);
  AA.resize(nSlv);
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
    if (mnpc[j] > 0)
    {
      AA(mnpc[j]) += elMat->b[0][j];
      gNA(mnpc[j]) += elMat->b[1][j];
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
  for (size_t j = 1; j <= AA.size(); j++)
    if (AA(j) > 0.0)
      gNA(j) /= AA(j);

#if INT_DEBUG > 2
  std::cout <<"MortarMats: Global AA:"<< AA
	    <<"MortarMats: Global gNA:"<< gNA
	    <<"MortarMats: Global phiA: "<< phiA;
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
  elm->A.front().resize(nedof,nen); // Element matrix phiA
  elm->b.front().resize(nen);       // Element vector AA
  elm->b.back().resize(nen);        // Element vector gNA
  return elm;
}


bool MortarContact::initElementBou (const std::vector<int>& MNPC,
				    LocalIntegral& elmInt)
{
  // Get nodal connectivities for the slave nodes of this element
  size_t nen = static_cast<ElmMats&>(elmInt).b.front().size();
  std::vector<int> MNPCu(MNPC.begin(),MNPC.begin()+nen);

  // Extract the current slave node displacements for this element
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
				      const Vector& lambda, double epsJxW) const
{
  size_t i, j, idof, jdof;
  usint  k, l;

  // Evaluate the geometric stiffness (Tg) of the contact spring
  // the normal vector (mNorm) of the master surface
  // and interpolation functions (Nm) of the master nodes
  Matrix Tg;
  master->geometricStiffness(X,mNorm,Tg,Nm);

  // Evaluate p_Nint at current integration point
  size_t nSlave = Ns.size();
  size_t nelnod = nSlave + Nm.size();
  double p_Nint = Ns.dot(lambda) * epsJxW;

  // Accumulate the displacement residuals

  for (j = jdof = 1; j <= nelnod; j++)
    for (l = 0; l < npv; l++, jdof++)
      R(jdof) += p_Nint*mNorm[l] * (j <= nSlave ? -Ns(j) : Nm(j-nSlave));

#if INT_DEBUG > 3
  std::cout <<"MortarContact::p_Nint = "<< p_Nint
	    <<"\nMortarContact::R"<< R;
#endif

  // Accumulate the tangent stiffness matrix (geometric stiffness terms only)

  if (Tg.normInf() <= 1.0e-12) return;

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

#if INT_DEBUG > 3
  std::cout <<"MortarContact::Kt"<< Kt;
#endif
}


bool MortarContact::assResAndTangent (SystemMatrix& Ktan, SystemVector& Res,
				      const MortarMats& mortar) const
{
  if (activeSlave.empty()) return true;

  std::vector<int> mnen, mnenj;
  size_t nnod = mortar.getNoNodes();
  size_t nslv = activeSlave.size();
  size_t i, j, ii, jj, n;
  usint  k, l;

  for (i = 1, ii = 0; i <= nnod; i++, ii += npv)
  {
    if (i <= mortar.getNoSlaves() && mortar.weightedArea(i) == 0.0)
      continue;
    else if (!mortar.getSAM().getNodeEqns(mnen,i))
      return false;

    // Add nodal sub-matrix on the diagonal
    Matrix eK(npv,npv);
    for (n = 1; n <= nslv; n++)
      if (activeSlave[n-1] && mortar.weightedArea(n) > 0.0)
      {
	double c = master->eps / mortar.weightedArea(n);
	for (k = 1; k <= npv; k++)
	  for (l = 1; l <= npv; l++)
	    eK(k,l) += c*mortar.phi(ii+k,n)*mortar.phi(ii+l,n);
      }

    if (!Ktan.assemble(eK,mortar.getSAM(),Res,mnen))
      return false;

    for (j = i+1, jj = ii+npv; j <= nnod; j++, jj += npv)
    {
      if (j <= mortar.getNoSlaves() && mortar.weightedArea(j) == 0.0)
	continue;
      else if (!mortar.getSAM().getNodeEqns(mnenj,j))
	return false;
      else if (mnen.size() == npv)
	mnen.insert(mnen.end(),mnenj.begin(),mnenj.end());
      else
	std::copy(mnenj.begin(),mnenj.end(),mnen.begin()+npv);

      // Add off-diagonal nodal sub-matrix symmetrically
      Matrix eM(npv+npv,npv+npv);
      for (n = 1; n <= nslv; n++)
	if (activeSlave[n-1] && mortar.weightedArea(n) > 0.0)
	{
	  double c = master->eps / mortar.weightedArea(n);
	  for (k = 1; k <= npv; k++)
	    for (l = 1; l <= npv; l++)
	      eM(k,npv+l) += c*mortar.phi(ii+k,n)*mortar.phi(jj+l,n);
	}

      // Lower triangular part
      for (k = 1; k <= npv; k++)
	for (l = 1; l <= npv; l++)
	  eM(npv+l,k) = eM(k,npv+l);

      if (!Ktan.assemble(eM,mortar.getSAM(),Res,mnen))
	return false;
    }
  }

  return true;
}


MortarPenalty::MortarPenalty (RigidBody* mst, const MortarMats& mats, usint nsd)
  : MortarContact(mst,NULL,nsd), mortar(mats)
{
  std::cout <<"\nMortarPenalty: Mortar contact formulation (Penalty, eps="
	    << master->eps <<")"<< std::endl;
}


void MortarPenalty::initIntegration (const TimeDomain& prm, const Vector&)
{
  // Initialize array of active slave nodes
  if (prm.first && prm.it == 0)
    activeSlave.resize(1,false); // Needed to initialize tangent matrix pattern
  else
    activeSlave.clear();

  activeSlave.reserve(mortar.getNoSlaves());

  size_t nmin = 0;
  double wmin = DBL_MAX;
  for (size_t n = 1; n <= mortar.getNoSlaves(); n++)
    if (mortar.weightedArea(n) > 0.0 && mortar.weightedGap(n) <= 0.0)
    {
      activeSlave.resize(n,false);
      activeSlave.back() = true;
    }
    else if (mortar.weightedArea(n) > 0.0 && mortar.weightedGap(n) < wmin)
    {
      nmin = n;
      wmin = mortar.weightedGap(n);
    }

#ifndef INT_DEBUG
  if (prm.it > 0) return;
#endif
  std::cout <<"  "<< std::count(activeSlave.begin(),activeSlave.end(),true)
	    <<" active contact nodes, closest non-active node "<< nmin
	    <<": wgap = "<< wmin << std::endl;
}


bool MortarPenalty::hasBoundaryTerms () const
{
  if (activeSlave.size() == 1)
    return activeSlave.front();

  return activeSlave.size() > 1;
}


LocalIntegral* MortarPenalty::getLocalIntegral (size_t nen, size_t, bool) const
{
  size_t nedof = npv*(nen+master->getNoNodes());
  ElmMats* elm = new ElmMats;
  elm->withLHS = true;
  elm->resize(1,2);
  elm->A.front().resize(nedof,nedof); // Element matrix Kt
  elm->b.front().resize(nedof);       // Element vector R
  elm->b.back().resize(nen);          // Element vector gNA
  return elm;
}


bool MortarPenalty::initElementBou (const std::vector<int>& MNPC,
				    LocalIntegral& elmInt)
{
  // Extract the current slave node displacements for this element
  if (!this->IntegrandBase::initElementBou(MNPC,elmInt))
    return false;

  // Extract the weighted nodal gaps for this element
  Vector& gNA = static_cast<ElmMats&>(elmInt).b.back();

  int lnod = 0;
  for (size_t j = 1; j <= gNA.size(); j++)
    if ((lnod = MNPC[j-1]) >= 0)
    {
      if ((size_t)lnod < nodMap.size())
	lnod = nodMap[lnod];
      else
	return false; // should not happen (logic error)

      gNA(j) = mortar.weightedGap(lnod);
      if (gNA(j) > 0.0) // we have separation at this point
	gNA(j) = 0.0;
    }

#if INT_DEBUG > 2
  std::cout <<"MortarPenalty::gNA"<< gNA;
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
			 elm.b.back(),master->eps*fe.detJxW);

  return true;
}


bool MortarPenalty::assemble (SystemMatrix& Ktan, SystemVector& Res) const
{
  size_t nnzK = Ktan.dim(0);
  bool locked = Ktan.lockPattern(false); // Unlock the sparsity pattern
  bool status = this->assResAndTangent(Ktan,Res,mortar);
  Ktan.lockPattern(locked); // Lock the sparsity pattern again
  size_t newnnz = Ktan.dim(0);
  if (newnnz > nnzK)
    std::cout <<"MortarPenalty: Tangent matrix has grown in size, from "
              << nnzK <<" to "<< newnnz <<" entries."<< std::endl;

  return status;
}


MortarAugmentedLag::MortarAugmentedLag (RigidBody* mst, const MortarMats& mats,
					const std::map<int,int>& alMap,
					usint nsd)
  : MortarContact(mst,NULL,nsd), mortar(mats), ALmap(alMap)
{
  std::cout <<"\nMortarAugmentedLag: Mortar contact formulation "
	    <<"(Augmented Lagrange, eps="<< master->eps <<")"<< std::endl;
}


void MortarAugmentedLag::initIntegration (const TimeDomain& prm,
                                          const Vector& psol)
{
  // Initialize array of active slave nodes and the Lagrange multipliers
  if (prm.first && prm.it == 0)
    activeSlave.resize(1,false); // Needed to initialize tangent matrix pattern
  else
    activeSlave.clear();

  activeSlave.reserve(mortar.getNoSlaves());
  lambda.resize(mortar.getNoSlaves(),true);

  const int* madof = mortar.getSAM().getMADOF();
  std::map<int,int>::const_iterator lit;

  size_t nmin = 0;
  double cmin = DBL_MAX;
  for (size_t n = 1; n <= lambda.size(); n++)
    if (mortar.weightedArea(n) > 0.0 && (lit = ALmap.find(n)) != ALmap.end())
    {
      // Get the Lagrange multiplier value at this slave node
      lambda(n) = psol(madof[lit->second-1]);
#if INT_DEBUG > 1
      std::cout <<"MortarAugmentedLag: Slave node "<< n
		<<" (ilag="<< madof[lit->second-1]
		<<"): wgap="<< mortar.weightedGap(n)
		<<", lambda="<< lambda(n) << std::endl;
#endif

      // Check the contact status at this slave node
      double cval = lambda(n) + master->eps*mortar.weightedGap(n);
      if (cval <= 0.0)
      {
	activeSlave.resize(n,false);
	activeSlave.back() = true;
      }
      else if (cval < cmin)
      {
	nmin = n;
	cmin = cval;
      }
    }

#ifndef INT_DEBUG
  if (prm.it > 0) return;
#endif
  std::cout <<"  "<< std::count(activeSlave.begin(),activeSlave.end(),true)
	    <<" active contact nodes, closest non-active node "<< nmin
	    <<": constraint = "<< cmin << std::endl;
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

  // Get nodal connectivities for the slave nodes of this element
  size_t nedof = elm.A.front().rows();
  size_t nen = nedof/npv - master->getNoNodes();
  std::vector<int> MNPCu(MNPC.begin(),MNPC.begin()+nen);

  // Extract the current slave node displacements for this element
  if (!this->IntegrandBase::initElementBou(MNPCu,elmInt))
    return false;

  // Extract the effective Lagrange multiplier values for this element
  Vector& LNA = elm.b.back();
  LNA.resize(nen);
  elm.iLag.resize(nen,0);

  int lnod = 0;
  size_t nLag = 0;
  for (size_t j = 1; j <= nen; j++)
    if ((lnod = MNPC[j-1]) >= 0)
    {
      nLag++;
      if ((size_t)lnod < nodMap.size())
	lnod = nodMap[lnod];
      else
	return false; // should not happen (logic error)

      if ((size_t)lnod <= activeSlave.size() && activeSlave[lnod-1])
      {
	elm.iLag[j-1] = nLag; // we have contact, store index for the Kul matrix
	LNA(j) = lambda(lnod) + master->eps*mortar.weightedGap(lnod);
      }
    }

  elm.b[0].resize(nedof+nLag); // Element vector R
  elm.A[1].resize(nedof,nLag); // Element matrix Kul

#if INT_DEBUG > 2
  std::cout <<"MortarAugmentedLag::LNA"<< LNA;
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
			 elm.b.back(),fe.detJxW);

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


bool MortarAugmentedLag::assemble (SystemMatrix& Ktan, SystemVector& Res) const
{
  size_t nnzK = Ktan.dim(0);
  bool locked = Ktan.lockPattern(false); // Unlock the sparsity pattern
  if (!this->assResAndTangent(Ktan,Res,mortar))
    return false;

  // Assemble the AL multiplier residual and associated tangent contribution
  Matrix  Kll(1,1);
  double& Kl = Kll(1,1);
  double  Rl, AA;

  double* Rvec = Res.getPtr();
  std::vector<int> mnen;

  // Loop over all Augmented Lagrange multipliers
  for (size_t n = 1; n <= mortar.getNoSlaves(); n++)
    if ((AA = mortar.weightedArea(n)) > 0.0)
    {
      if (n <= activeSlave.size() && activeSlave[n-1])
      {
	// Residual contribution for the active AL multiplier
	Rl = -mortar.weightedGap(n) * AA;
	Kl = 0.0;
      }
      else
      {
	// Residual and tangent contribution for the inactive AL multiplier
	Rl = (lambda(n)/master->eps) * AA;
	Kl = master->eps * AA;
      }

      // Find the global equation number of this AL multiplier
      std::map<int,int>::const_iterator lit = ALmap.find(n);
      if (lit == ALmap.end())
         return false; // logic error
      else if (!mortar.getSAM().getNodeEqns(mnen,lit->second))
        return false;

#if INT_DEBUG > 3
      std::cout <<"MortarAugmentedLag::assemble: ieq="<< mnen.front()
		<<": Rl="<< Rl <<" Kll="<< Kl << std::endl;
#endif

      // Add into the global residual vector tangent stiffness matrix
      Rvec[mnen.front()-1] += Rl;
      if (Kl != 0.0)
	Ktan.assemble(Kll,mortar.getSAM(),Res,mnen);
    }

  Res.restore(Rvec);
  Ktan.lockPattern(locked); // Lock the sparsity pattern again
  size_t newnnz = Ktan.dim(0);
  if (newnnz >  nnzK)
    std::cout <<"MortarAugmentedLag: Tangent matrix has grown in size, from "
              << nnzK <<" to "<< newnnz <<" entries."<< std::endl;

  return true;
}


MortarAugmentedLag::ALElmMats::ALElmMats (size_t nedof)
{
  withLHS = true;
  this->resize(3,2);
  A.front().resize(nedof,nedof);
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

    // Insert the coupling matrix Kul symmetrically
    for (j = 1; j <= nlag; j++)
    {
      K(i,ndof+j) = A[1](i,j);
      K(ndof+j,i) = A[1](i,j);
    }
  }

#if INT_DEBUG > 2
  std::cout <<"\nMortarAugmentedLag::ALElmMats: Newton matrix ="<< A.back();
#endif
  return A.back();
}
