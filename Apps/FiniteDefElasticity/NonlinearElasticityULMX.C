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
#include "Tensor.h"
#include "Vec3Oper.h"
#include "int_debug.h"

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
  void acckmx2d_(const int& axS, const int& nEN, const double* Nr,
		 const double* dNdx, const double* dNdxBar,
		 const double* D, double* eKt);
  //! \brief Accumulates material stiffness contributions for mixed 3D problems.
  void acckmx3d_(const int& nEN, const double* dNdx, const double* dNdxBar,
		 const double* D, double* eKt);
}
#endif


/*!
  \brief A struct with integration point data needed by \a finalizeElement.
*/

struct ItgPtData
{
  Tensor F;     //!< Deformation gradient at current step/iteration
  Tensor Fp;    //!< Deformation gradient at previous (converged) step
  Vector Nr;    //!< Basis function values (for axisymmetric problems)
  Matrix dNdx;  //!< Basis function gradients at current configuration
  Vector Phi;   //!< Internal modes for the pressure/volumetric-change fields
  Vec4   X;     //!< Cartesian coordinates of current integration point
  double detJW; //!< Jacobian determinant times integration point weight
  //! \brief Default constructor.
  ItgPtData() : F(3), Fp(3) { detJW = 0.0; }
};


/*!
  \brief Class containing internal data for a mixed element.
*/

class MxElmData
{
public:
  //! \brief Default constructor.
  MxElmData() : Hh(0), iP(0) {}
  //! \brief Empty destructor.
  virtual ~MxElmData() {}

  //! \brief Initializes the element data.
  bool init(const Vec3& Xc, size_t nPt, int p, unsigned short int nsd)
  {
    iP = 0;
    X0 = Xc;
    myData.resize(nPt);
    if (Hh)
    {
      size_t nPm = utl::Pascal(p,nsd);
      Hh->resize(nPm,nPm,true);
      return true;
    }
    std::cerr <<" *** MxElmData::init: No Hh matrix defined."<< std::endl;
    return false;
  }

  Vec3    X0; //!< Cartesian coordinates of the element center
  Matrix* Hh; //!< Pointer to element Hh-matrix associated with pressure modes
  size_t  iP; //!< Local integration point counter

  std::vector<ItgPtData> myData; //!< Local integration point data
};


/*!
  \brief Class collecting element matrices associated with a mixed FEM problem.
*/

class MxMats : public MxElmData, public ElmMats {};


/*!
  \brief Class representing integrated norm quantities for a mixed element.
*/

class MxNorm : public MxElmData, public LocalIntegral
{
public:
  //! \brief The constructor initializes the ElmNorm pointer.
  MxNorm(LocalIntegral& n)
  {
    myNorm = static_cast<ElmNorm*>(&n);
    Hh = new Matrix();
  }

  //! \brief Alternative constructor allocating an internal ElmNorm.
  MxNorm(size_t n)
  {
    myNorm = new ElmNorm(n);
    Hh = new Matrix();
  }

  //! \brief The destructor destroys the ElmNorm object.
  virtual ~MxNorm() { myNorm->destruct(); }

  //! \brief Returns the LocalIntegral object to assemble into the global one.
  virtual const LocalIntegral* ref() const { return myNorm; }

  ElmNorm* myNorm; //!< Pointer to the actual element norms object
};


NonlinearElasticityULMX::NonlinearElasticityULMX (unsigned short int n,
						  bool axS, int pp)
  : NonlinearElasticityUL(n,axS), p(pp)
{
  // Both the current and previous solutions are needed
  primsol.resize(2);
}


void NonlinearElasticityULMX::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityULMX: Mixed formulation,"
	    <<" discontinuous pressure, p="<< p << std::endl;

  this->NonlinearElasticityUL::print(os);
}


LocalIntegral* NonlinearElasticityULMX::getLocalIntegral (size_t nen, size_t,
							  bool neumann) const
{
  MxMats* result = new MxMats;

  switch (m_mode)
    {
    case SIM::STATIC:
      result->withLHS = !neumann;
    case SIM::RHS_ONLY:
      result->rhsOnly = neumann;
      result->resize(neumann?1:2,1);
      break;

    case SIM::DYNAMIC:
      result->withLHS = !neumann;
      result->rhsOnly = neumann;
      result->resize(neumann?1:3,1);
      break;

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->resize(1,0);
      break;

    default:
      return result;
    }

  result->Hh = &result->A.back();

  for (size_t i = 1; i < result->A.size(); i++)
    result->A[i-1].resize(nsd*nen,nsd*nen);

  if (result->b.size())
    result->b.front().resize(nsd*nen);

  return result;
}


bool NonlinearElasticityULMX::initElement (const std::vector<int>& MNPC,
					   const Vec3& Xc, size_t nPt,
                                           LocalIntegral& elmInt)
{
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityULMX::initElement";
  std::cout <<", Xc = "<< Xc << std::endl;
#endif

  return static_cast<MxMats&>(elmInt).init(Xc,nPt,p,nsd) &&
         this->IntegrandBase::initElement(MNPC,elmInt);
}


bool NonlinearElasticityULMX::evalInt (LocalIntegral& elmInt,
				       const FiniteElement& fe,
				       const TimeDomain&, const Vec3& X) const
{
  MxElmData* mx = 0;
  MxNorm* mxNrm = 0;
  MxMats* mxMat = dynamic_cast<MxMats*>(&elmInt);
  if (mxMat)
    mx = mxMat;
  else if ((mxNrm = dynamic_cast<MxNorm*>(&elmInt)))
    mx = mxNrm;
  else
    return false;

  ItgPtData& ptData = mx->myData[mx->iP++];
  const Vectors& eV = mxNrm ? mxNrm->myNorm->vec : elmInt.vec;
  if (eV.size() < 2) return false;


#if INT_DEBUG > 1
  std::cout <<"NonlinearElasticityULMX::dNdX ="<< fe.dNdX;
#endif

  // Evaluate the deformation gradient, Fp, at previous configuration
  if (!this->formDefGradient(eV[1],fe.N,fe.dNdX,X.x,ptData.Fp))
    return false;

  // Evaluate the deformation gradient, F, at current configuration
  if (!this->formDefGradient(eV[0],fe.N,fe.dNdX,X.x,ptData.F))
    return false;

  if (nDF == 2) // In 2D we always assume plane strain so set F(3,3)=1
    ptData.F(3,3) = ptData.Fp(3,3) = 1.0;

  // Invert the deformation gradient ==> Fi
  Matrix Fi(nsd,nsd);
  if (nsd == 3)
    Fi.fill(ptData.F.ptr());
  else
    for (unsigned short int i = 1; i <= nsd; i++)
      for (unsigned short int j = 1; j <= nsd; j++)
	Fi(i,j) = ptData.F(i,j);

  double J = Fi.inverse();
  if (axiSymmetry) J *= ptData.F(3,3);
  if (J == 0.0) return false;

  // Push-forward the basis function gradients to current configuration
  ptData.dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
  if (axiSymmetry)
    ptData.Nr = fe.N * (1.0/(X.x + eV.front().dot(fe.N,0,nsd)));

  ptData.X.assign(X);
  ptData.detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Evaluate the pressure modes (generalized coordinates)
  Vec3 Xg = X - mx->X0;
  if (nsd == 3)
    utl::Pascal(p,Xg.x,Xg.y,Xg.z,ptData.Phi);
  else if (nsd == 2)
    utl::Pascal(p,Xg.x,Xg.y,ptData.Phi);
  else
    return false;

  if (mx->Hh)
    // Integrate the Hh-matrix
    mx->Hh->outer_product(ptData.Phi,ptData.Phi*ptData.detJW,true);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(mxMat->A[eM-1],fe.N,X,J*ptData.detJW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(mxMat->b[eS-1],fe.N,X,J*ptData.detJW);

  return true;
}


bool NonlinearElasticityULMX::finalizeElement (LocalIntegral& elmInt,
					       const TimeDomain& prm, size_t iG)
{
  if (!iS && !eKm && !eKg) return true;

  MxMats& mx = static_cast<MxMats&>(elmInt);
  if (!mx.Hh) return false;

  size_t iP;
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityULMX::finalizeElement\n";
  for (iP = 1; iP <= mx.myData.size(); iP++)
  {
    const ItgPtData& pt = mx.myData[iP-1];
    std::cout <<"\n     X"   << iP <<" = "<< pt.X;
    std::cout <<"\n     detJ"<< iP <<"W = "<< pt.detJW;
    std::cout <<"\n     F"   << iP <<" =\n"<< pt.F;
    std::cout <<"     dNdx"<< iP <<" ="<< pt.dNdx;
    std::cout <<"     Phi" << iP <<" ="<< pt.Phi;
  }
  std::cout <<"\n     H ="<< *mx.Hh << std::endl;
#endif

  // 1. Eliminate the internal pressure DOFs by static condensation.

  size_t i, j;
  size_t nPM = mx.Hh->rows();
  size_t nEN = mx.myData.front().dNdx.rows();
  size_t nGP = mx.myData.size();

  // Invert the H-matrix
  if (!utl::invert(*mx.Hh)) return false;

  const Matrix& Hi = *mx.Hh;

  Matrix Theta(2,nGP), dNdxBar[nGP];

  if (nPM == 1)
  {
    // Simplified calculation when the only pressure mode is Phi=1.0

    double h1 = 0.0;
    double h2 = 0.0;
    dNdxBar->resize(nEN,nsd,true);

    for (iP = 0; iP < nGP; iP++)
    {
      const ItgPtData& pt = mx.myData[iP];
      double dVol = pt.F.det()*pt.detJW;
      h1 += dVol;
      h2 += pt.Fp.det()*pt.detJW;
      dNdxBar->add(pt.dNdx,dVol);
      if (axiSymmetry)
	for (j = 1; j <= nEN; j++)
	  dNdxBar[0](j,1) += pt.Nr(j)*dVol;
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
      const ItgPtData& pt = mx.myData[iP];
      for (i = 0; i < nPM; i++)
      {
	double h1 = pt.Phi[i] * pt.F.det() * pt.detJW;
	double h2 = pt.Phi[i] * pt.Fp.det()* pt.detJW;
	J1[i] += h1;
	J2[i] += h2;
	Gg[i].add(pt.dNdx,h1);
	if (axiSymmetry)
	  for (j = 1; j <= nEN; j++)
	    Gg[i](j,1) += pt.Nr(j)*h1;
      }
    }

#if INT_DEBUG > 1
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

#if INT_DEBUG > 1
    std::cout <<"NonlinearElasticityULMX::Hj1 ="<< Hj1;
    std::cout <<"NonlinearElasticityULMX::Hj2 ="<< Hj2;
    for (i = 0; i < nPM; i++)
      std::cout <<"NonlinearElasticityULMX::Hg"<< i+1 <<" ="<< Hg[i];
#endif

    for (iP = 0; iP < nGP; iP++)
    {
      i = iP+1;
      ItgPtData& pt = mx.myData[iP];

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
    ItgPtData& pt = mx.myData[iP-1];
    detF(iP) = pt.F.det();
    pt.F  *= pow(fabs(Theta(1,iP)/detF(iP)),1.0/3.0);
    pt.Fp *= pow(fabs(Theta(2,iP)/pt.Fp.det()),1.0/3.0);
  }

#if INT_DEBUG > 1
  std::cout <<"NonlinearElasticityULMX::detF ="<< detF;
  std::cout <<"NonlinearElasticityULMX::Theta ="<< Theta;
  for (iP = 0; iP < nGP; iP++)
    std::cout <<"NonlinearElasticityULMX::dNdxBar"<< iP+1 <<" ="<< dNdxBar[iP];
#endif

  // 2. Evaluate the constitutive relation and calculate mixed pressure state.

  Matrix Cmat, D[nGP];
  Vector Sigm(nPM), Hsig;
  std::vector<SymmTensor> Sig;
  Sig.reserve(nGP);

  bool lHaveStress = false;
  for (iP = 0; iP < nGP; iP++)
  {
#if INT_DEBUG > 0
    std::cout <<"\n   iGP =     "<< iP+1 << std::endl;
#endif
    const ItgPtData& pt = mx.myData[iP];

    // Evaluate the constitutive relation
    double U = 0.0;
    Sig.push_back(SymmTensor(3));
    if (!material->evaluate(Cmat,Sig.back(),U,iG++,pt.X,pt.F,
			    SymmTensor(0),1,&prm))
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
    if (!mx.Hh->multiply(Sigm,Hsig)) return false;

  // 3. Evaluate tangent matrices.

  double Press = 0.0, Bpres = 0.0;
  SymmTensor Sigma(nsd,axiSymmetry);
  for (iP = 0; iP < nGP; iP++)
  {
    const ItgPtData& pt = mx.myData[iP];
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
	this->formKG(mx.A[eKg-1],pt.Nr,pt.dNdx,1.0,Sigma,dFmix*pt.detJW);
    }

    if (eKm)
    {
#ifdef USE_FTNMAT
      mdma3d_(Bpres,Press,Sig[iP].ptr(),D[iP].ptr(),INT_DEBUG,6);

      // Integrate the material stiffness matrix
      D[iP] *= dFmix*pt.detJW;
      if (nsd == 2)
	acckmx2d_(axiSymmetry,nEN,pt.Nr.ptr(),pt.dNdx.ptr(),dNdxBar[iP].ptr(),
		  D[iP].ptr(),mx.A[eKm-1].ptr());
      else
	acckmx3d_(nEN,pt.dNdx.ptr(),dNdxBar[iP].ptr(),
		  D[iP].ptr(),mx.A[eKm-1].ptr());
#endif
    }

    if (iS && lHaveStress)
    {
      // Compute the small-deformation strain-displacement matrix B from dNdx
      Matrix Bmat;
      if (axiSymmetry)
	this->formBmatrix(Bmat,pt.Nr,pt.dNdx,1.0);
      else
	this->formBmatrix(Bmat,pt.dNdx);

      // Integrate the internal forces
      Sigma *= -dFmix*pt.detJW;
      if (!Bmat.multiply(Sigma,mx.b[iS-1],true,true)) // ES -= B^T*sigma
	return false;
    }
  }

#if INT_DEBUG > 0
  std::cout <<"\n *** Leaving NonlinearElasticityULMX::finalizeElement";
  if (eKm) std::cout <<", eKm ="<< mx.A[eKm-1];
  if (iS) std::cout <<" iS ="<< mx.b[iS-1];
  std::cout << std::endl;
#endif
  return true;
}


NormBase* NonlinearElasticityULMX::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormULMX(*const_cast<NonlinearElasticityULMX*>(this));
}


LocalIntegral* ElasticityNormULMX::getLocalIntegral (size_t nen, size_t iEl,
						     bool neumann) const
{
  LocalIntegral* result = this->NormBase::getLocalIntegral(nen,iEl,neumann);
  if (result) return new MxNorm(static_cast<ElmNorm&>(*result));

  // Element norms are not requested, so allocate an inter object instead that
  // will delete itself when invoking the destruct method.
  return new MxNorm(this->getNoFields());
}


bool ElasticityNormULMX::initElement (const std::vector<int>& MNPC,
				      const Vec3& Xc, size_t nPt,
				      LocalIntegral& elmInt)
{
  NonlinearElasticityULMX& p = static_cast<NonlinearElasticityULMX&>(myProblem);
  MxNorm& mx = static_cast<MxNorm&>(elmInt);

  return mx.init(Xc,nPt,p.p,p.nsd) &&
         this->NormBase::initElement(MNPC,*mx.myNorm);
}


bool ElasticityNormULMX::initElementBou (const std::vector<int>& MNPC,
					 LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,*static_cast<MxNorm&>(elmInt).myNorm);
}


bool ElasticityNormULMX::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const TimeDomain& td, const Vec3& X) const
{
  NonlinearElasticityULMX& p = static_cast<NonlinearElasticityULMX&>(myProblem);

  return p.evalInt(elmInt,fe,td,X);
}


bool ElasticityNormULMX::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X, const Vec3& normal) const
{
  return this->ElasticityNormUL::evalBou(*static_cast<MxNorm&>(elmInt).myNorm,
					 fe,X,normal);
}


bool ElasticityNormULMX::finalizeElement (LocalIntegral& elmInt,
					  const TimeDomain& prm, size_t iG)
{
  NonlinearElasticityULMX& p = static_cast<NonlinearElasticityULMX&>(myProblem);
  MxNorm& mx = static_cast<MxNorm&>(elmInt);

  if (!mx.Hh) return false;

  size_t iP;
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering ElasticityNormULMX::finalizeElement\n";
  for (iP = 1; iP <= mx.myData.size(); iP++)
  {
    const ItgPtData& pt = mx.myData[iP-1];
    std::cout <<"\n     X"   << iP <<" = "<< pt.X;
    std::cout <<"\n     detJ"<< iP <<"W = "<< pt.detJW;
    std::cout <<"\n     F"   << iP <<" =\n"<< pt.F;
    std::cout <<"     dNdx"<< iP <<" ="<< pt.dNdx;
    std::cout <<"     Phi" << iP <<" ="<< pt.Phi;
  }
  std::cout <<"\n     H ="<< *(mx.Hh) << std::endl;
#endif

  // 1. Eliminate the internal pressure DOFs by static condensation.

  size_t nPM = mx.Hh->rows();
  size_t nGP = mx.myData.size();

  // Invert the H-matrix
  Matrix& Hi = *mx.Hh;
  if (!utl::invert(Hi)) return false;

  Vector Theta(nGP);

  if (nPM == 1)
  {
    // Simplified calculation when the only pressure mode is Phi=1.0

    double h1 = 0.0;
    for (iP = 0; iP < nGP; iP++)
      h1 += mx.myData[iP].F.det() * mx.myData[iP].detJW;

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
      const ItgPtData& pt = mx.myData[iP];
      for (size_t i = 0; i < nPM; i++)
	Ji[i] += pt.Phi[i] * pt.F.det() * pt.detJW;
    }

    // Calculate Hj = Hi*Ji
    if (!Hi.multiply(Ji,Hj)) return false;

#if INT_DEBUG > 1
    std::cout <<"ElasticityNormULMX::Ji ="<< Ji;
    std::cout <<"ElasticityNormULMX::Hj ="<< Hj;
#endif

    // Calculate Theta = Hj*Phi
    for (iP = 0; iP < nGP; iP++)
      Theta[iP] = Hj.dot(mx.myData[iP].Phi);
  }

#if INT_DEBUG > 1
  std::cout <<"ElasticityNormULMX::Theta ="<< Theta;
#endif

  // 2. Evaluate the constitutive relation and integrate energy.

  Matrix Cmat;
  SymmTensor Sig(3), Eps(3);
  for (iP = 0; iP < nGP; iP++)
  {
#if INT_DEBUG > 0
    std::cout <<"\n   iGP =     "<< iP+1 << std::endl;
#endif
    ItgPtData& pt = mx.myData[iP];
    Eps.rightCauchyGreen(pt.F); // Green Lagrange strain tensor
    Eps -= 1.0;
    Eps *= 0.5;

    // Modify the deformation gradient
    Tensor Fbar(pt.F);
    Fbar *= pow(fabs(Theta[iP]/pt.F.det()),1.0/3.0);

    // Compute the strain energy density, U(Eps) = Int_Eps (S:E) dE
    // and the Cauchy stress tensor, Sig
    double U = 0.0;
    if (!p.material->evaluate(Cmat,Sig,U,iG++,pt.X,Fbar,Eps,3,&prm,&pt.F))
      return false;

    // Integrate the norms
    if (!ElasticityNormUL::evalInt(*mx.myNorm,Sig,U,Theta[iP],pt.detJW))
      return false;
  }

  return true;
}
