// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityULMixed.C
//!
//! \date Dec 22 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity mixed problems.
//!
//==============================================================================

#include "NonlinearElasticityULMixed.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ElmMats.h"
#include "ElmNorm.h"
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
  //! \brief Accumulates material stiffness contributions for 2D problems.
  void acckm2d_(const int& nEN, const double* dNdx,
		const double* D, double* eKt);
  //! \brief Accumulates material stiffness contributions for 3D problems.
  void acckm3d_(const int& nEN, const double* dNdx,
		const double* D, double* eKt);
}
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif


//! \brief Enum for element-level solution vectors.
enum SolutionVectors {
  U = 0, // displacement
  T = 1, // volumetric change (theta)
  P = 2, // pressure
  NSOL = 3
};

//! \brief Enum for element-level right-hand-side vectors.
enum ResidualVectors {
  Ru = 0,
  Rt = 1,
  Rp = 2,
  Rres = 3,
  NVEC = 4
};

//! \brief Enum for element-level left-hand-side matrices.
enum TangentMatrices {
  Kuu = 0,
  Kut = 1,
  Kup = 2,
  Ktt = 3,
  Ktp = 4,
  Ktan = 5,
  NMAT = 6
};


NonlinearElasticityULMixed::MixedElmMats::MixedElmMats ()
{
  this->resize(NMAT,NVEC);
}


/*!
  Tangent matrix for the mixed formulation:
  \f[{\bf K}_{\rm T} = \left[\begin{array}{ccc}
  {\bf K}_{uu} & {\bf K}_{u\theta} & {\bf K}_{up} \\
  {\bf K}_{u\theta}^T & {\bf K}_{\theta\theta} & {\bf K}_{\theta p} \\
  {\bf K}_{up}^T & {\bf K}_{\theta p}^T & {\bf 0}
  \end{array}\right]\f]
  where
  -# \f${\bf K}_{uu}\f$ :
     tangent stiffness matrix of the displacement formulation
  -# \f${\bf K}_{\theta\theta}\f$ :
     tangent stiffness matrix related to volumetric change
  -# \f${\bf K}_{u\theta}\f$ :
     coupling stiffness between displacements and volumetric change
  -# \f${\bf K}_{up}\f$ :
     coupling stiffness between displacements and pressure
  -# \f${\bf K}_{\theta p}\f$ :
     coupling stiffness between volumetric change and pressure
*/

const Matrix& NonlinearElasticityULMixed::MixedElmMats::getNewtonMatrix () const
{
#if INT_DEBUG > 0
  std::cout <<"\nNonlinearElasticityULMixed::MixedElmMats::getNewtonMatrix:"
	    <<"\nKuu ="<< A[Kuu]
	    <<"\nKut ="<< A[Kut]
	    <<"\nKup ="<< A[Kup]
	    <<"\nKtt ="<< A[Ktt]
	    <<"\nKtp ="<< A[Ktp];
#endif

  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  size_t i, j;
  size_t n = A[Kuu].rows();
  size_t m = A[Ktt].rows();

  for (i = 1; i <= n; i++)
  {
    for (j = 1; j <= n; j++)
      N(i,j) = A[Kuu](i,j);
    for (j = 1; j <= m; j++)
    {
      size_t k = n+2*j-1;
      size_t l = n+2*j;
      N(i,k) = A[Kut](i,j);
      N(k,i) = A[Kut](i,j);
      N(i,l) = A[Kup](i,j);
      N(l,i) = A[Kup](i,j);
    }
  }

  for (i = 1; i <= m; i++)
    for (j = 1; j <= m; j++)
    {
      size_t ki = n+2*i-1;
      size_t kj = n+2*j-1;
      size_t lj = n+2*j;
      N(ki,kj) = A[Ktt](i,j);
      N(ki,lj) = A[Ktp](i,j);
      N(lj,ki) = A[Ktp](i,j);
    }

#if INT_DEBUG > 0
  std::cout <<"\nNewton matrix ="<< A[Ktan];
#endif
  return A[Ktan];
}


/*!
  Right-hand-side vector for the mixed formulation:
  \f[{\bf R} = \left\{\begin{array}{c}
  {\bf R}_u \\ {\bf R}_\theta \\ {\bf R}_p
  \end{array}\right\}\f]
  where
  -# \f${\bf R}_{u}\f$ :
     residual forces of the displacement DOFs
  -# \f${\bf R}_\theta\f$ :
     residual forces of the volumetric change DOFs
  -# \f${\bf R}_p\f$ :
     residual forces of the pressure DOFs
*/

const Vector& NonlinearElasticityULMixed::MixedElmMats::getRHSVector () const
{
#if INT_DEBUG > 0
  std::cout <<"\nNonlinearElasticityULMixed::MixedElmMats::getRHSVector:"
	    <<"\nRu ="<< b[Ru];
  if (withLHS)
    std::cout <<"\nRt ="<< b[Rt] <<"\nRp ="<< b[Rp];
#endif
  Vector& R = const_cast<Vector&>(b[Rres]);

  size_t i;
  size_t n = b[Ru].size();
  size_t m = b[Rt].size();

  for (i = 1; i <= n; i++)
    R(i) = b[Ru](i);

  if (withLHS)
    for (i = 1; i <= m; i++)
    {
      R(n+2*i-1) = b[Rt](i);
      R(n+2*i  ) = b[Rp](i);
    }
  else
    std::fill(R.begin()+n,R.end(),0.0);

#if INT_DEBUG > 0
  std::cout <<"\nRHS ="<< b[Rres];
#endif
  return b[Rres];
}


NonlinearElasticityULMixed::NonlinearElasticityULMixed (unsigned short int n)
  : NonlinearElasticityUL(n), Fbar(3), Dmat(7,7)
{
  if (myMats) delete myMats;
  myMats = new MixedElmMats();
  mySols.resize(NSOL);

  eKm = &myMats->A[Kuu];
  eKg = &myMats->A[Kuu];

  iS = &myMats->b[Ru];
  eS = &myMats->b[Ru];
  eV = &mySols[U];
}


void NonlinearElasticityULMixed::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityULMixed: "
	    <<"Continuous volumetric change and pressure fields"<< std::endl;

  this->NonlinearElasticityUL::print(os);
}


void NonlinearElasticityULMixed::setMode (SIM::SolutionMode mode)
{
  switch (mode)
    {
    case SIM::INIT:
    case SIM::STATIC:
      tracVal.clear();
      myMats->rhsOnly = false;
      break;

    case SIM::RHS_ONLY:
      tracVal.clear();
    case SIM::RECOVERY:
      myMats->rhsOnly = true;
      break;

    default:
      std::cerr <<"\n *** NonlinearElasticityULMixed::setMode: Invalid mode "
		<< mode << std::endl;
    }
}


static int iP = 0; //!< Local integration point counter for debug output
std::vector<int> mixedDbgEl; //!< List of elements for additional output


bool NonlinearElasticityULMixed::initElement (const std::vector<int>& MNPC1,
					      const std::vector<int>& MNPC2,
					      size_t n1)
{
  iP = 0;

  const size_t nen1 = MNPC1.size();
  const size_t nen2 = MNPC2.size();
  const size_t nedof = nsd*nen1 + 2*nen2;

  myMats->A[Kut].resize(nsd*nen1,nen2,true);
  myMats->A[Kup].resize(nsd*nen1,nen2,true);
  myMats->A[Ktt].resize(nen2,nen2,true);
  myMats->A[Ktp].resize(nen2,nen2,true);
  myMats->A[Ktan].resize(nedof,nedof,true);

  myMats->b[Rt].resize(nen2,true);
  myMats->b[Rp].resize(nen2,true);
  myMats->b[Rres].resize(nedof,true);

  // Extract the element level volumetric change and pressure vectors
  int ierr = 0;
  if (!primsol.front().empty())
    if ((ierr = utl::gather(MNPC2,0,2,primsol.front(),mySols[T],nsd*n1,n1) +
		utl::gather(MNPC2,1,2,primsol.front(),mySols[P],nsd*n1,n1)))
      std::cerr <<" *** NonlinearElasticityULMixed::initElement: Detected "
		<< ierr/2 <<" node numbers out of range."<< std::endl;

  // The other element matrices are initialized by the parent class method
  return this->NonlinearElasticityUL::initElement(MNPC1) && ierr == 0;
}


bool NonlinearElasticityULMixed::initElementBou (const std::vector<int>& MNPC1,
						 const std::vector<int>& MNPC2,
						 size_t)
{
  const size_t nen1 = MNPC1.size();
  const size_t nen2 = MNPC2.size();
  const size_t nedof = nsd*nen1 + 2*nen2;

  myMats->b[Ru].resize(nsd*nen1,true);
  myMats->b[Rres].resize(nedof,true);

  // Extract the element level displacement vector
  int ierr = 0;
  if (!primsol.front().empty())
    if ((ierr = utl::gather(MNPC1,nsd,primsol.front(),mySols[U])))
      std::cerr <<" *** NonlinearElasticityULMixed::initElementBou: Detected "
                << ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


bool NonlinearElasticityULMixed::evalIntMx (LocalIntegral*& elmInt,
					    const MxFiniteElement& fe,
					    const TimeDomain& prm,
					    const Vec3& X) const
{
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityULMixed::evalIntMx: iP = "
	    << ++iP <<", (X,t) = "<< X << std::endl;
#endif

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(fe.dN1dX,F,E))
    return false;

  bool lHaveStrains = !E.isZero(1.0e-16);

  // Evaluate the volumetric change and pressure fields
  double Theta = fe.N2.dot(mySols[T]) + 1.0;
  double Press = fe.N2.dot(mySols[P]);
#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityULMixed::b_theta ="<< mySols[T];
  std::cout <<"NonlinearElasticityULMixed::b_p ="<< mySols[P];
  std::cout <<"NonlinearElasticityULMixed::Theta = "<< Theta
	    <<", Press = "<< Press << std::endl;
#endif

  // Compute the mixed model deformation gradient, F_bar
  Fbar = F; // notice that F_bar always has dimension 3
  double r1 = pow(fabs(Theta/F.det()),1.0/3.0);

  Fbar *= r1;
  if (nsd == 2) // In 2D we always assume plane strain so set F(3,3)=1
    Fbar(3,3) = r1;
  else if (nsd != 3)
    return false;

  // Invert the deformation gradient ==> Fi
  Matrix Fi(nsd,nsd);
  Fi.fill(F.ptr());
  double J = Fi.inverse();
  if (J == 0.0) return false;

  // Push-forward the basis function gradients to current configuration
  dNdx.multiply(fe.dN1dX,Fi); // dNdx = dNdX * F^-1

  // Compute the mixed integration point volume
  double dVol = Theta*fe.detJxW;
  double dVup = J*fe.detJxW;

#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityULMixed::dNdX ="<< fe.dN1dX;
  std::cout <<"NonlinearElasticityULMixed::dNdx ="<< dNdx;
  std::cout <<"NonlinearElasticityULMixed::Fbar =\n"<< Fbar;
  std::cout <<"NonlinearElasticityULMixed::E =\n"<< E;
  std::cout <<"NonlinearElasticityULMixed::dVol = "<< dVol
	    <<" "<< dVup << std::endl;
#endif

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,fe.N1,X,J*fe.detJxW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,fe.N1,X,J*fe.detJxW);

  // Evaluate the constitutive relation
  SymmTensor Sig(3), Sigma(nsd);
  double U, Bpres = 0.0, Mpres = 0.0;
  if (!material->evaluate(Cmat,Sig,U,X,Fbar,E,lHaveStrains,&prm))
    return false;

#ifdef USE_FTNMAT
  // Project the constitutive matrix for the mixed model
  pcst3d_(Cmat.ptr(),Dmat.ptr(),INT_DEBUG,6);
#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityULMixed::Dmat ="<< Dmat;
#endif
#else
  return false;
#endif

  if (lHaveStrains)
  {
    if (prm.it == 0 &&
	find(mixedDbgEl.begin(),mixedDbgEl.end(),fe.iel) != mixedDbgEl.end())
    {
#if INT_DEBUG > 0
      if (iP == 1)
#else
      if (++iP == 1)
#endif
	std::cout <<"\n  ** Stresses in element "<< fe.iel << std::endl;

      std::cout << iP;
      const RealArray& sigma = Sig;
      for (size_t i = 0; i < sigma.size(); i++)
	std::cout <<" "<< sigma[i];
      std::cout <<" "<< Sig.vonMises() << std::endl;
    }

    // Compute mixed strees
    Bpres = Sig.trace()/3.0;
    Mpres = Press * J/Theta;
    Sigma = Sig + (Mpres-Bpres);

    // Integrate the geometric stiffness matrix
    this->formKG(*eKg,dNdx,Sigma,dVol);
  }

#ifdef USE_FTNMAT
  mdma3d_(Bpres,Mpres,Sig.ptr(),Dmat.ptr(),INT_DEBUG,6);

  // Integrate the material stiffness matrix
  Dmat *= dVol;
  if (nsd == 2)
    acckm2d_(fe.N1.size(),dNdx.ptr(),Dmat.ptr(),eKm->ptr());
  else
    acckm3d_(fe.N1.size(),dNdx.ptr(),Dmat.ptr(),eKm->ptr());
#endif

  // Integrate the volumetric change and pressure tangent terms
  size_t a, b;
  unsigned short int i, j, k;
  for (i = 1; i <= 6; i++)
    Dmat(i,7) /= Theta;
  Dmat(7,7) /= (Theta*Theta);
  for (a = 1; a <= fe.N1.size(); a++)
    for (b = 1; b <= fe.N2.size(); b++)
      for (i = 1; i <= nsd; i++)
      {
	myMats->A[Kut](nsd*(a-1)+i,b) += dNdx(a,i)*fe.N2(b) * Dmat(i,7);
	if (nsd == 2)
	  myMats->A[Kut](nsd*(a-1)+i,b) += dNdx(a,3-i)*fe.N2(b) * Dmat(4,7);
	else if (i < 3)
	{
	  j = i + 1;
	  k = j%3 + 1;
	  myMats->A[Kut](nsd*(a-1)+i,b) += dNdx(a,3-i)*fe.N2(b) * Dmat(4,7);
	  myMats->A[Kut](nsd*(a-1)+j,b) += dNdx(a,5-j)*fe.N2(b) * Dmat(5,7);
	  myMats->A[Kut](nsd*(a-1)+k,b) += dNdx(a,4-k)*fe.N2(b) * Dmat(6,7);
	}
	myMats->A[Kup](nsd*(a-1)+i,b) += dNdx(a,i)*fe.N2(b) * dVup;
      }

  myMats->A[Ktt].outer_product(fe.N2,fe.N2*Dmat(7,7),true);    // += N2*N2^T*D77
  myMats->A[Ktp].outer_product(fe.N2,fe.N2*(-fe.detJxW),true); // -= N2*N2^T*|J|

  if (lHaveStrains)
  {
    // Compute the small-deformation strain-displacement matrix B from dNdx
    if (!this->Elasticity::formBmatrix(dNdx))
      return false;

    // Integrate the internal forces
    Sigma *= -dVol;
    if (!Bmat.multiply(Sigma,myMats->b[Ru],true,true))  // -= B^T*sigma
      return false;

    // Integrate the volumetric change and pressure forces
    myMats->b[Rt].add(fe.N2,(Press - Bpres)*fe.detJxW); // += N2*(p-pBar)*|J|
    myMats->b[Rp].add(fe.N2,(Theta - J)*fe.detJxW);     // += N2*(Theta-J)*|J|

#if INT_DEBUG > 4
    std::cout <<"NonlinearElasticityULMixed::Sigma*dVol =\n"<< Sigma;
    std::cout <<"NonlinearElasticityULMixed::Bmat ="<< Bmat;
    std::cout <<"NonlinearElasticityULMixed::Ru("<< iP <<") ="<< myMats->b[Ru];
#endif
  }

  return this->getIntegralResult(elmInt);
}


/*!
  The boundary integral is the same as that of the parent class. It does not
  depend on the pressure and volumetric change fields. Thus, this call is
  forwarded to the corresponding single-field method of the parent class.
*/

bool NonlinearElasticityULMixed::evalBouMx (LocalIntegral*& elmInt,
					    const MxFiniteElement& fe,
					    const Vec3& X,
					    const Vec3& normal) const
{
  return this->evalBou(elmInt,fe,X,normal);
}


size_t NonlinearElasticityULMixed::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd + 2;
  else
    return this->Elasticity::getNoFields(fld);
}


const char* NonlinearElasticityULMixed::getField1Name (size_t i,
						       const char* prefix) const
{
  if (i < nsd || (i > (size_t)(nsd+1) && i != 12))
    return this->Elasticity::getField1Name(i,prefix);

  static std::string name;
  if (prefix)
    name = prefix + std::string(" ");
  else
    name.clear();

  if (i == 12)
    name += "theta+p";
  else if (i == nsd)
    name += "theta";
  else
    name += "p";

  return name.c_str();
}


NormBase* NonlinearElasticityULMixed::getNormIntegrand (AnaSol*) const
{
  NonlinearElasticityULMixed* ulp;
  ulp = const_cast<NonlinearElasticityULMixed*>(this);
  return new ElasticityNormULMixed(*ulp);
}


bool ElasticityNormULMixed::evalIntMx (LocalIntegral*& elmInt,
				       const MxFiniteElement& fe,
				       const TimeDomain& prm,
				       const Vec3& X) const
{
  ElmNorm& pnorm = NormBase::getElmNormBuffer(elmInt);

  NonlinearElasticityULMixed* ulp;
  ulp = static_cast<NonlinearElasticityULMixed*>(&myProblem);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Tensor F(fe.dN1dX.cols());
  SymmTensor E(F.dim());
  if (!ulp->kinematics(fe.dN1dX,F,E))
    return false;

  // Evaluate the volumetric change field
  double Theta = fe.N2.dot(ulp->mySols[T]) + 1.0;

  // Compute the mixed model deformation gradient, F_bar
  ulp->Fbar = F; // notice that F_bar always has dimension 3
  double r1 = pow(fabs(Theta/F.det()),1.0/3.0);

  ulp->Fbar *= r1;
  if (F.dim() == 2) // In 2D we always assume plane strain so set F(3,3)=1
    ulp->Fbar(3,3) = r1;

  // Compute the strain energy density, U(E) = Int_E (Sig:Eps) dEps
  double U = 0.0;
  SymmTensor Sig(3);
  if (!ulp->material->evaluate(ulp->Cmat,Sig,U,X,ulp->Fbar,E,3,&prm,&F))
    return false;

  // Integrate the energy norm a(u^h,u^h) = Int_Omega0 U(E) dV0
  pnorm[0] += U*fe.detJxW;
  // Integrate the L2-norm ||Sig|| = Int_Omega0 Sig:Sig dV0
  pnorm[2] += Sig.L2norm(false)*fe.detJxW;
  // Integrate the von Mises stress norm
  pnorm[3] += Sig.vonMises(false)*fe.detJxW;

  return true;
}


bool ElasticityNormULMixed::evalBouMx (LocalIntegral*& elmInt,
				       const MxFiniteElement& fe,
				       const Vec3& X, const Vec3& normal) const
{
  return this->ElasticityNormUL::evalBou(elmInt,fe,X,normal);
}
