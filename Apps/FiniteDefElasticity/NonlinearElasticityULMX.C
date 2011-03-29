// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityULMX.C
//!
//! \date Dec 14 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity mixed problems.
//!
//==============================================================================

#include "NonlinearElasticityULMX.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3Oper.h"

#ifdef USE_FTNMAT
extern "C" {
  //! \brief Project the constitutive matrix for mixed models.
  void pcst3d_(const double* C, double* D,
	       const int& ipsw, const int& iwr);
  //! \brief Compute finite deformation mixed tangent material matrix.
  void mdma3d_(const double& P_bar, const double& P_mix,
	       const double* Sig, double* D,
	       const int& ipsw, const int& iwr);
  //! \brief Accumulates material stiffness contributions for mixed 2D problems.
  void acckmx2d_(const int& nEN, const double* dNdx, const double* dNdxBar,
		 const double* D, double* eKt);
  //! \brief Accumulates material stiffness contributions for mixed 3D problems.
  void acckmx3d_(const int& nEN, const double* dNdx, const double* dNdxBar,
		 const double* D, double* eKt);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


NonlinearElasticityULMX::NonlinearElasticityULMX (unsigned short int n, int pp)
  : NonlinearElasticityUL(n)
{
  p = pp;
  iP = 0;

  // Both the current and previous solutions are needed
  primsol.resize(2);
}


void NonlinearElasticityULMX::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityULMX: Mixed formulation,"
	    <<" discontinuous pressure, p="<< p << std::endl;

  this->NonlinearElasticityUL::print(os);
}


void NonlinearElasticityULMX::setMode (SIM::SolutionMode mode)
{
  if (!myMats) return;

  myMats->rhsOnly = false;
  Hh = eM = eKm = eKg = 0;
  iS = eS = eV  = 0;

  switch (mode)
    {
    case SIM::STATIC:
      myMats->resize(2,3);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      Hh  = &myMats->A[1];
      break;

    case SIM::DYNAMIC:
      myMats->resize(3,3);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      eM  = &myMats->A[1];
      Hh  = &myMats->A[2];
      break;

    case SIM::RHS_ONLY:
      if (myMats->A.size() < 2)
	myMats->resize(2,3);
      else
	myMats->b.resize(3);
      eKm = &myMats->A[0];
      eKg = &myMats->A[0];
      Hh  = &myMats->A.back();
      myMats->rhsOnly = true;
      break;

    case SIM::RECOVERY:
      if (myMats->A.size() < 1)
	myMats->resize(1,3);
      else
	myMats->b.resize(3);
      Hh = &myMats->A.back();
      eV = &myMats->b[1];
      myMats->rhsOnly = true;
      return;

    default:
      this->Elasticity::setMode(mode);
      return;
    }

  // We always need the force+displacement vectors in nonlinear simulations
  iS = &myMats->b[0];
  eS = &myMats->b[0];
  eV = &myMats->b[1];
  tracVal.clear();
}


bool NonlinearElasticityULMX::initElement (const std::vector<int>& MNPC,
					   const Vec3& Xc, size_t nPt)
{
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityULMX::initElement";
  std::cout <<", Xc = "<< Xc << std::endl;
#endif

  iP = 0;
  X0 = Xc;
  myData.resize(nPt);
  size_t nPm = utl::Pascal(p,nsd);
  if (Hh) Hh->resize(nPm,nPm,true);

  // The other element matrices are initialized by the parent class method
  return this->NonlinearElasticityUL::initElement(MNPC);
}


bool NonlinearElasticityULMX::evalInt (LocalIntegral*& elmInt,
				       const FiniteElement& fe,
				       const TimeDomain&, const Vec3& X) const
{
  if (myMats->b.size() < 3) return false;

  ItgPtData& ptData = myData[iP++];

#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityUL::dNdX ="<< fe.dNdX;
#endif

  // Evaluate the deformation gradient, Fp, at previous configuration
  const_cast<NonlinearElasticityULMX*>(this)->eV = &myMats->b[2];
  if (!this->formDefGradient(fe.dNdX,ptData.Fp))
    return false;

  // Evaluate the deformation gradient, F, at current configuration
  const_cast<NonlinearElasticityULMX*>(this)->eV = &myMats->b[1];
  if (!this->formDefGradient(fe.dNdX,ptData.F))
    return false;

  if (nsd == 2) // In 2D we always assume plane strain so set F(3,3)=1
    ptData.F(3,3) = ptData.Fp(3,3) = 1.0;

  // Invert the deformation gradient ==> Fi
  Matrix Fi(nsd,nsd);
  for (unsigned short int i = 1; i <= nsd; i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      Fi(i,j) = ptData.F(i,j);
  double J = Fi.inverse();
  if (J == 0.0) return false;

  // Push-forward the basis function gradients to current configuration
  ptData.dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1

  ptData.X.assign(X);
  ptData.detJW = fe.detJxW;

  // Evaluate the pressure modes (generalized coordinates)
  Vec3 Xg = X-X0;
  if (nsd == 3)
    utl::Pascal(p,Xg.x,Xg.y,Xg.z,ptData.Phi);
  else if (nsd == 2)
    utl::Pascal(p,Xg.x,Xg.y,ptData.Phi);
  else
    return false;

  if (Hh)
    // Integrate the Hh-matrix
    Hh->outer_product(ptData.Phi,ptData.Phi*fe.detJxW,true);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,fe.N,X,J*fe.detJxW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,fe.N,X,J*fe.detJxW);

  return this->getIntegralResult(elmInt);
}


bool NonlinearElasticityULMX::finalizeElement (LocalIntegral*& elmInt,
					       const TimeDomain& prm)
{
  if (!iS && !eKm && !eKg) return true;
  if (!Hh) return false;

#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityULMX::finalizeElement\n";
  for (iP = 1; iP <= myData.size(); iP++)
  {
    const ItgPtData& pt = myData[iP-1];
    std::cout <<"\n     X"   << iP <<" = "<< pt.X;
    std::cout <<"\n     detJ"<< iP <<"W = "<< pt.detJW;
    std::cout <<"\n     F"   << iP <<" =\n"<< pt.F;
    std::cout <<"     dNdx"<< iP <<" ="<< pt.dNdx;
    std::cout <<"     Phi" << iP <<" ="<< pt.Phi;
  }
  std::cout <<"\n     H ="<< *Hh << std::endl;
#endif

  // 1. Eliminate the internal pressure DOFs by static condensation.

  size_t i, j;
  size_t nPM = Hh->rows();
  size_t nEN = myData.front().dNdx.rows();
  size_t nGP = myData.size();

  // Invert the H-matrix
  if (!utl::invert(*Hh)) return false;

  const Matrix& Hi = *Hh;

  Matrix Theta(2,nGP), dNdxBar[nGP];

  if (nPM == 1)
  {
    // Simplified calculation when the only pressure mode is Phi=1.0

    double h1 = 0.0;
    double h2 = 0.0;
    dNdxBar->resize(nEN,nsd,true);

    for (iP = 0; iP < nGP; iP++)
    {
      const ItgPtData& pt = myData[iP];
      h1 += pt.F.det() * pt.detJW;
      h2 += pt.Fp.det()* pt.detJW;
      dNdxBar->add(pt.dNdx,pt.F.det()*pt.detJW);
    }
    dNdxBar->multiply(1.0/h1);

    // All gauss point values are identical
    Theta(1,1) = h1 * Hi(1,1);
    Theta(2,1) = h2 * Hi(1,1);
    for (iP = 2; iP <= nGP; iP++)
    {
      Theta(1,iP) = Theta(1,1);
      Theta(2,iP) = Theta(2,1);
      dNdxBar[iP-1] = dNdxBar[0];
    }
  }
  else // higher order elements with nPM > 1
  {
    Matrix Gg[nPM], Hg[nPM];
    Vector J1(nPM), J2(nPM), Hj1, Hj2;
    for (i = 0; i < nPM; i++)
    {
      Gg[i].resize(nEN,nsd,true);
      Hg[i].resize(nEN,nsd,true);
    }

    // Integrate the Ji- and G-matrices
    for (iP = 0; iP < nGP; iP++)
    {
      const ItgPtData& pt = myData[iP];
      for (i = 0; i < nPM; i++)
      {
	double h1 = pt.Phi[i] * pt.F.det() * pt.detJW;
	double h2 = pt.Phi[i] * pt.Fp.det()* pt.detJW;
	J1[i] += h1;
	J2[i] += h2;
	Gg[i].add(pt.dNdx,h1);
      }
    }

#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityULMX::J1 ="<< J1;
    std::cout <<"NonlinearElasticityULMX::J2 ="<< J2;
    for (i = 0; i < nPM; i++)
      std::cout <<"NonlinearElasticityULMX::G"<< i+1 <<" ="<< Gg[i];
#endif

    // Calculate Hji = Hi*Ji
    if (!Hi.multiply(J1,Hj1)) return false;
    if (!Hi.multiply(J2,Hj2)) return false;

    // Calculate Hg = Hi*Gg
    for (i = 1; i <= nPM; i++)
      for (j = 1; j <= nPM; j++)
	Hg[i-1].add(Gg[j-1],Hi(i,j));

#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityULMX::Hj1 ="<< Hj1;
    std::cout <<"NonlinearElasticityULMX::Hj2 ="<< Hj2;
    for (i = 0; i < nPM; i++)
      std::cout <<"NonlinearElasticityULMX::Hg"<< i+1 <<" ="<< Hg[i];
#endif

    for (iP = 0; iP < nGP; iP++)
    {
      i = iP+1;
      ItgPtData& pt = myData[iP];

      // Calculate Theta = Hj*Phi
      Theta(1,i) = Hj1.dot(pt.Phi);
      Theta(2,i) = Hj2.dot(pt.Phi);

      // Calculate dNdxBar = Hg*Phi * 1/Theta
      dNdxBar[iP].resize(nEN,nsd,true);
      for (j = 0; j < nPM; j++)
	dNdxBar[iP].add(Hg[j],pt.Phi[j]/Theta(1,i));
    }
  }

  // Modify the deformation gradients
  Vector detF(nGP);
  for (iP = 1; iP <= nGP; iP++)
  {
    ItgPtData& pt = myData[iP-1];
    detF(iP) = pt.F.det();
    pt.F  *= pow(fabs(Theta(1,iP)/detF(iP)),1.0/3.0);
    pt.Fp *= pow(fabs(Theta(2,iP)/pt.Fp.det()),1.0/3.0);
  }

#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityULMX::detF ="<< detF;
  std::cout <<"NonlinearElasticityULMX::Theta ="<< Theta;
  for (iP = 0; iP < nGP; iP++)
    std::cout <<"NonlinearElasticityULMX::dNdxBar"<< iP+1 <<" ="<< dNdxBar[iP];
#endif

  // 2. Evaluate the constitutive relation and calculate mixed pressure state.

  Matrix D[nGP];
  Vector Sigm(nPM), Hsig;
  std::vector<SymmTensor> Sig;
  Sig.reserve(nGP);

  bool lHaveStress = false;
  for (iP = 0; iP < nGP; iP++)
  {
#if INT_DEBUG > 0
    std::cout <<"\n   iGP =     "<< iP+1 << std::endl;
#endif
    const ItgPtData& pt = myData[iP];

    // Evaluate the constitutive relation
    double U = 0.0;
    Sig.push_back(SymmTensor(3));
    if (!material->evaluate(Cmat,Sig.back(),U,pt.X,pt.F,Sig.back(),true,&prm))
      return false;

#ifdef USE_FTNMAT
    // Project the constitutive matrix for the mixed model
    D[iP].resize(7,7);
    pcst3d_(Cmat.ptr(),D[iP].ptr(),INT_DEBUG,6);
#else
    std::cerr <<" *** NonlinearElasticityULMX::finalizeElement: "
	      <<" Does not work when compiled without USE_FTNMAT"<< std::endl;
    return false;
#endif

    if (!Sig.back().isZero())
    {
      // Integrate mean pressure
      lHaveStress = true;
      double Mpress = Sig.back().trace()*pt.detJW/3.0;
      for (i = 1; i <= nPM; i++)
	Sigm(i) += Mpress * pt.Phi(i);
    }
  }

  // Divide pressure by reference volume; Hsig = Hi*Sigm
  if (lHaveStress)
    if (!Hh->multiply(Sigm,Hsig)) return false;

  // 3. Evaluate tangent matrices.

  double Press = 0.0, Bpres = 0.0;
  SymmTensor Sigma(nsd);
  for (iP = 0; iP < nGP; iP++)
  {
    const ItgPtData& pt = myData[iP];
    double dFmix = Theta(1,iP+1);
    if (lHaveStress)
    {
      Press = Hsig.dot(pt.Phi) * detF[iP]/dFmix;
      Bpres = Sig[iP].trace()/3.0;
      Sigma = Sig[iP] + (Press-Bpres);
#if INT_DEBUG > 0
      std::cout <<"NonlinearElasticityULMX::Sigma"<< iP+1 <<"\n"<< Sigma;
#endif
      if (eKg)
	// Integrate the geometric stiffness matrix
	this->formKG(*eKg,pt.dNdx,Sigma,dFmix*pt.detJW);
    }

    if (eKm)
    {
#ifdef USE_FTNMAT
      mdma3d_(Bpres,Press,Sig[iP].ptr(),D[iP].ptr(),INT_DEBUG,6);

      // Integrate the material stiffness matrix
      D[iP] *= dFmix*pt.detJW;
      if (nsd == 2)
	acckmx2d_(nEN,pt.dNdx.ptr(),dNdxBar[iP].ptr(),D[iP].ptr(),eKm->ptr());
      else
	acckmx3d_(nEN,pt.dNdx.ptr(),dNdxBar[iP].ptr(),D[iP].ptr(),eKm->ptr());
#endif
    }

    if (iS && lHaveStress)
    {
      // Compute the small-deformation strain-displacement matrix B from dNdx
      if (!this->formBmatrix(pt.dNdx)) return false;

      // Integrate the internal forces
      Sigma *= -dFmix*pt.detJW;
      if (!Bmat.multiply(Sigma,*iS,true,true)) // ES -= B^T*sigma
	return false;
    }
  }

#if INT_DEBUG > 0
  std::cout <<"\n *** Leaving NonlinearElasticityULMX::finalizeElement";
  if (eKm) std::cout <<", eKm ="<< *eKm;
  if (iS) std::cout <<" iS ="<< *iS;
  std::cout << std::endl;
#endif
  return true;
}


NormBase* NonlinearElasticityULMX::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormULMX(*const_cast<NonlinearElasticityULMX*>(this));
}


bool ElasticityNormULMX::initElement (const std::vector<int>& MNPC,
				      const Vec3& Xc, size_t nPt)
{
  NonlinearElasticityULMX& elp = static_cast<NonlinearElasticityULMX&>(problem);

  elp.iP = 0;
  elp.X0 = Xc;
  elp.myData.resize(nPt);
  size_t nPm = utl::Pascal(elp.p,elp.nsd);
  if (elp.Hh) elp.Hh->resize(nPm,nPm,true);

  // The other element matrices are initialized by the parent class method
  return elp.NonlinearElasticityUL::initElement(MNPC);
}


bool ElasticityNormULMX::evalInt (LocalIntegral*& elmInt,
				  const FiniteElement& fe,
				  const TimeDomain& td, const Vec3& X) const
{
  return static_cast<NonlinearElasticityULMX&>(problem).evalInt(elmInt,fe,td,X);
}


bool ElasticityNormULMX::finalizeElement (LocalIntegral*& elmInt,
					  const TimeDomain& prm)
{
  NonlinearElasticityULMX& elp = static_cast<NonlinearElasticityULMX&>(problem);
  if (!elp.Hh) return false;

  size_t iP;
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering ElasticityNormULMX::finalizeElement\n";
  for (iP = 1; iP <= elp.myData.size(); iP++)
  {
    const NonlinearElasticityULMX::ItgPtData& pt = elp.myData[iP-1];
    std::cout <<"\n     X"   << iP <<" = "<< pt.X;
    std::cout <<"\n     detJ"<< iP <<"W = "<< pt.detJW;
    std::cout <<"\n     F"   << iP <<" =\n"<< pt.F;
    std::cout <<"     dNdx"<< iP <<" ="<< pt.dNdx;
    std::cout <<"     Phi" << iP <<" ="<< pt.Phi;
  }
  std::cout <<"\n     H ="<< *(elp.Hh) << std::endl;
#endif

  // 1. Eliminate the internal pressure DOFs by static condensation.

  size_t nPM = elp.Hh->rows();
  size_t nGP = elp.myData.size();

  // Invert the H-matrix
  Matrix& Hi = *elp.Hh;
  if (!utl::invert(Hi)) return false;

  Vector Theta(nGP);

  if (nPM == 1)
  {
    // Simplified calculation when the only pressure mode is Phi=1.0

    double h1 = 0.0;
    for (iP = 0; iP < nGP; iP++)
      h1 += elp.myData[iP].F.det() * elp.myData[iP].detJW;

    // All gauss point values are identical
    Theta.front() = h1 * Hi(1,1);
    for (iP = 1; iP < nGP; iP++)
      Theta[iP] = Theta.front();
  }
  else // higher order elements with nPM > 1
  {
    Vector Ji(nPM), Hj;

    // Integrate the Ji-matrix
    for (iP = 0; iP < nGP; iP++)
    {
      const NonlinearElasticityULMX::ItgPtData& pt = elp.myData[iP];
      for (size_t i = 0; i < nPM; i++)
	Ji[i] += pt.Phi[i] * pt.F.det() * pt.detJW;
    }

    // Calculate Hj = Hi*Ji
    if (!Hi.multiply(Ji,Hj)) return false;

#if INT_DEBUG > 0
    std::cout <<"ElasticityNormULMX::Ji ="<< Ji;
    std::cout <<"ElasticityNormULMX::Hj ="<< Hj;
#endif

    // Calculate Theta = Hj*Phi
    for (iP = 0; iP < nGP; iP++)
      Theta[iP] = Hj.dot(elp.myData[iP].Phi);
  }

#if INT_DEBUG > 0
  std::cout <<"ElasticityNormULMX::Theta ="<< Theta;
#endif

  // 2. Evaluate the constitutive relation and integrate energy.

  ElmNorm& pnorm = ElasticityNorm::getElmNormBuffer(elmInt);

  Matrix C;
  SymmTensor Sig(3);

  for (iP = 0; iP < nGP; iP++)
  {
#if INT_DEBUG > 0
    std::cout <<"\n   iGP =     "<< iP+1 << std::endl;
#endif
    NonlinearElasticityULMX::ItgPtData& pt = elp.myData[iP];

    // Modify the deformation gradient
    pt.F *= pow(fabs(Theta[iP]/pt.F.det()),1.0/3.0);

    // Compute the strain energy density
    double U = 0.0;
    if (!elp.material->evaluate(C,Sig,U,pt.X,pt.F,Sig,false,&prm))
      return false;

    // Integrate strain energy
    pnorm[0] += U*pt.detJW;
  }

  return true;
}
