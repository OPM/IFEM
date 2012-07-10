// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityFbar.C
//!
//! \date Aug 18 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity with Fbar method.
//!
//==============================================================================

#include "NonlinearElasticityFbar.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "GaussQuadrature.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "Tensor.h"
#include "int_debug.h"


/*!
  \brief A struct with volumetric sampling point data.
*/

struct VolPtData
{
  double J;    //!< Determinant of current deformation gradient
  Vector Nr;   //!< Basis function values (for axisymmetric problems)
  Matrix dNdx; //!< Basis function gradients at current configuration
  //! \brief Default constructor.
  VolPtData() { J = 1.0; }
};


/*!
  \brief Class containing internal data for an Fbar element.
*/

class FbarElmData
{
public:
  //! \brief Default constructor.
  FbarElmData() : pbar(0), scale(0.0), iP(0) {}
  //! \brief Empty destructor.
  virtual ~FbarElmData() {}

  //! \brief Initializes the element data
  bool init(unsigned short int nsd, size_t nPt)
  {
    pbar = ceil(pow((double)nPt,1.0/(double)nsd)-0.5);
    if (pow((double)pbar,(double)nsd) == (double)nPt)
      myVolData.resize(nPt);
    else
    {
      pbar = 0;
      std::cerr <<" *** FbarElmData::init: Invalid element, "<< nPt
		<<" volumetric sampling points specified."<< std::endl;
      return false;
    }
    scale = pbar > 1 ? 1.0/GaussQuadrature::getCoord(pbar)[pbar-1] : 1.0;
    iP = 0;
    return true;
  }

  int    pbar;  //!< Polynomial order of the internal volumetric data field
  double scale; //!< Scaling factor for extrapolation from sampling points
  size_t iP;    //!< Volumetric sampling point counter

  std::vector<VolPtData> myVolData; //!< Volumetric sampling point data
};


/*!
  \brief Class collecting element matrices associated with a Fbar FEM problem.
*/

class FbarMats : public FbarElmData, public ElmMats {};


/*!
  \brief Class representing integrated norm quantities for an Fbar element.
*/

class FbarNorm : public FbarElmData, public LocalIntegral
{
public:
  //! \brief The constructor initializes the ElmNorm pointer.
  FbarNorm(LocalIntegral& n) { myNorm = static_cast<ElmNorm*>(&n); }
  //! \brief Alternative constructor allocating an internal ElmNorm.
  FbarNorm(size_t n) { myNorm = new ElmNorm(n); }
  //! \brief The destructor destroys the ElmNorm object.
  virtual ~FbarNorm() { myNorm->destruct(); }

  //! \brief Returns the LocalIntegral object to assemble into the global one.
  virtual const LocalIntegral* ref() const { return myNorm; }

  ElmNorm* myNorm; //!< Pointer to the actual element norms object
};


NonlinearElasticityFbar::NonlinearElasticityFbar (unsigned short int n,
						  bool axS, int nvp)
  : NonlinearElasticityUL(n,axS), npt1(nvp)
{
}


void NonlinearElasticityFbar::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityFbar: F-bar formulation, "<< npt1
	    <<" volumetric points in each direction."<< std::endl;

  this->NonlinearElasticityUL::print(os);
}


LocalIntegral* NonlinearElasticityFbar::getLocalIntegral (size_t nen, size_t,
							  bool neumann) const
{
  ElmMats* result = new FbarMats;
  switch (m_mode)
  {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:1,1);
      break;

    case SIM::DYNAMIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:2,1);
      break;

    case SIM::VIBRATION:
    case SIM::BUCKLING:
      result->withLHS = true;
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
    case SIM::MASS_ONLY:
      result->withLHS = true;
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->rhsOnly = true;
      result->resize(neumann?0:1,1);
      break;

    case SIM::RECOVERY:
      result->rhsOnly = true;
      break;

    default:
      ;
  }

  for (size_t i = 0; i < result->A.size(); i++)
    result->A[i].resize(nsd*nen,nsd*nen);

  if (result->b.size())
    result->b.front().resize(nsd*nen);

  return result;
}


bool NonlinearElasticityFbar::initElement (const std::vector<int>& MNPC,
                                           const Vec3&, size_t nPt,
					   LocalIntegral& elmInt)
{
  FbarMats& fbar = static_cast<FbarMats&>(elmInt);

#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityFbar::initElement";
  std::cout <<", nPt = "<< nPt <<", pbar = "<< fbar.pbar
	    <<", scale = "<< fbar.scale << std::endl;
#endif
  return fbar.init(nsd,nPt) && this->IntegrandBase::initElement(MNPC,elmInt);
}


bool NonlinearElasticityFbar::reducedInt (LocalIntegral& elmInt,
					  const FiniteElement& fe,
					  const Vec3& X) const
{
  FbarNorm* fbNorm = 0;
  FbarElmData* fbar = dynamic_cast<FbarMats*>(&elmInt);
  if (!fbar)
  {
    if ((fbNorm = dynamic_cast<FbarNorm*>(&elmInt)))
      fbar = fbNorm;
    else
      return false;
  }

#if INT_DEBUG > 1
  std::cout <<"NonlinearElasticityFbar::u(red) = "<< fe.u;
  if (nsd > 1) std::cout <<" "<< fe.v;
  if (nsd > 2) std::cout <<" "<< fe.w;
  std::cout <<" (iP="<< fbar->iP+1 <<")\n"
	    <<"NonlinearElasticityFbar::dNdX ="<< fe.dNdX;
#endif

  const Vector& eV = fbNorm ? fbNorm->myNorm->vec.front() : elmInt.vec.front();
  VolPtData& ptData = fbar->myVolData[fbar->iP++];

  // Evaluate the deformation gradient, F, at current configuration
  Matrix B;
  Tensor F(nDF);
  SymmTensor E(nsd,axiSymmetry);
  if (!this->kinematics(eV,fe.N,fe.dNdX,X.x,B,F,E))
    return false;

  if (E.isZero(1.0e-16))
  {
    // Initial state, no deformation yet
    ptData.J = 1.0;
    ptData.dNdx = fe.dNdX;
    if (axiSymmetry && X.x > 0.0)
      ptData.Nr = fe.N * (1.0/X.x);
  }
  else
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    if (nsd == nDF)
      Fi.fill(F.ptr());
    else
      for (unsigned short int i = 1; i <= nsd; i++)
	for (unsigned short int j = 1; j <= nsd; j++)
	  Fi(i,j) = F(i,j);

    ptData.J = Fi.inverse();
    if (axiSymmetry) ptData.J *= F(3,3);
    if (ptData.J == 0.0)
      return false;

    // Push-forward the basis function gradients to current configuration
    ptData.dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
    if (axiSymmetry && X.x > 0.0)
      ptData.Nr = fe.N * (1.0/(X.x + eV.dot(fe.N,0,nsd)));

#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityFbar::J = "<< ptData.J
	      <<"\nNonlinearElasticityFbar::dNdx ="<< ptData.dNdx;
    if (axiSymmetry)
      std::cout <<"NonlinearElasticityFbar::Nr ="<< ptData.Nr;
#endif
  }
#if INT_DEBUG > 0
  if (fbar->iP == fbar->myVolData.size())
    std::cout <<"NonlinearElasticityFbar: Volumetric sampling points completed."
	      <<"\n"<< std::endl;
#endif

  return true;
}


/*!
  \brief Computes the discrete gradient operator, \b G, for continuum elements.
*/

static void getGradOperator (double r, const Vector& N, const Matrix& dNdx,
			     Matrix& G)
{
  size_t nen = dNdx.rows();
  size_t nsd = dNdx.cols();
  size_t nrg = r > 0.0 ? nsd*nsd+1 : nsd*nsd;
  G.resize(nrg,nsd*nen);

  size_t a, i, j, k, l = 0;
  for (a = 1; a <= nen; a++, l += nsd)
  {
    for (i = k = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++, k++)
	G(k,l+j) = dNdx(a,i);
    if (k == nrg)
      G(k,l+1) = N(a)/r;
  }
}


/*!
  \brief Computes the spatial tangent modulus matrix, \b A.
  \details For 2D continuum elements.
*/

static void getAmat2D (const Matrix& C, const SymmTensor& Sig, Matrix& A)
{
  A.resize(4,4);

  A(1,1) = C(1,1) + Sig(1,1);
  A(2,1) = C(3,1);
  A(3,1) = C(3,1) + Sig(2,1);
  A(4,1) = C(2,1);

  A(1,2) = C(1,3);
  A(2,2) = C(3,3) + Sig(1,1);
  A(3,2) = C(3,3);
  A(4,2) = C(2,3) + Sig(1,2);

  A(1,3) = C(1,3) + Sig(1,2);
  A(2,3) = C(3,3);
  A(3,3) = C(3,3) + Sig(2,2);
  A(4,3) = C(2,3);

  A(1,4) = C(1,2);
  A(2,4) = C(3,2) + Sig(2,1);
  A(3,4) = C(3,2);
  A(4,4) = C(2,2) + Sig(2,2);
}


/*!
  \brief Computes the spatial tangent modulus matrix, \b A.
  \details For 3D axisymmetric continuum elements.
*/

static void getAmatAx (const Matrix& C, const SymmTensor& Sig, Matrix& A)
{
  A.resize(5,5);

  A(1,1) = C(1,1) + Sig(1,1);
  A(2,1) = C(4,1);
  A(3,1) = C(4,1) + Sig(2,1);
  A(4,1) = C(2,1);
  A(5,1) = C(3,1);

  A(1,2) = C(1,4);
  A(2,2) = C(4,4) + Sig(1,1);
  A(3,2) = C(4,4);
  A(4,2) = C(2,4) + Sig(1,2);
  A(5,2) = C(3,4);

  A(1,3) = C(1,4) + Sig(1,2);
  A(2,3) = C(4,4);
  A(3,3) = C(4,4) + Sig(2,2);
  A(4,3) = C(2,4);
  A(5,3) = C(3,4);

  A(1,4) = C(1,2);
  A(2,4) = C(4,2) + Sig(2,1);
  A(3,4) = C(4,2);
  A(4,4) = C(2,2) + Sig(2,2);
  A(5,4) = C(3,2);

  A(1,5) = C(1,3);
  A(2,5) = C(4,3);
  A(3,5) = C(4,3);
  A(4,5) = C(2,3);
  A(5,5) = C(3,3) + Sig(3,3);
}


/*!
  \brief Computes the spatial tangent modulus matrix, \b A.
  \details For 3D continuum elements.
*/

static void getAmat3D (const Matrix& C, const SymmTensor& Sig, Matrix& A)
{
  A.resize(9,9);

  A(1,1) = C(1,1) + Sig(1,1);
  A(2,1) = C(4,1);
  A(3,1) = C(6,1);
  A(4,1) = C(4,1) + Sig(2,1);
  A(5,1) = C(2,1);
  A(6,1) = C(5,1);
  A(7,1) = C(6,1) + Sig(3,1);
  A(8,1) = C(5,1);
  A(9,1) = C(3,1);

  A(1,2) = C(1,4);
  A(2,2) = C(4,4) + Sig(1,1);
  A(3,2) = C(6,4);
  A(4,2) = C(4,4);
  A(5,2) = C(2,4) + Sig(2,1);
  A(6,2) = C(5,4);
  A(7,2) = C(6,4);
  A(8,2) = C(5,4) + Sig(3,1);
  A(9,2) = C(3,4);

  A(1,3) = C(1,6);
  A(2,3) = C(4,6);
  A(3,3) = C(6,6) + Sig(1,1);
  A(4,3) = C(4,6);
  A(5,3) = C(2,6);
  A(6,3) = C(5,6) + Sig(1,2);
  A(7,3) = C(6,6);
  A(8,3) = C(5,6);
  A(9,3) = C(3,6) + Sig(1,2);

  A(1,4) = C(1,4) + Sig(1,2);
  A(2,4) = C(4,4);
  A(3,4) = C(6,4);
  A(4,4) = C(4,4) + Sig(2,2);
  A(5,4) = C(2,4);
  A(6,4) = C(5,4);
  A(7,4) = C(6,4) + Sig(3,2);
  A(8,4) = C(5,4);
  A(9,4) = C(3,4);

  A(1,5) = C(1,2);
  A(2,5) = C(4,2) + Sig(1,2);
  A(3,5) = C(6,2);
  A(4,5) = C(4,2);
  A(5,5) = C(2,2) + Sig(2,2);
  A(6,5) = C(5,2);
  A(7,5) = C(6,2);
  A(8,5) = C(5,2) + Sig(3,2);
  A(9,5) = C(3,2);

  A(1,6) = C(1,5);
  A(2,6) = C(4,5);
  A(3,6) = C(6,5) + Sig(2,1);
  A(4,6) = C(4,5);
  A(5,6) = C(2,5);
  A(6,6) = C(5,5) + Sig(2,2);
  A(7,6) = C(6,5);
  A(8,6) = C(5,5);
  A(9,6) = C(3,5) + Sig(2,3);

  A(1,7) = C(1,6) + Sig(1,3);
  A(2,7) = C(4,6);
  A(3,7) = C(6,6);
  A(4,7) = C(4,6) + Sig(2,3);
  A(5,7) = C(2,6);
  A(6,7) = C(5,6);
  A(7,7) = C(6,6) + Sig(3,3);
  A(8,7) = C(5,6);
  A(9,7) = C(3,6);

  A(1,8) = C(1,5);
  A(2,8) = C(4,5) + Sig(1,3);
  A(3,8) = C(6,5);
  A(4,8) = C(4,5);
  A(5,8) = C(2,5) + Sig(2,3);
  A(6,8) = C(5,5);
  A(7,8) = C(6,5);
  A(8,8) = C(5,5) + Sig(3,3);
  A(9,8) = C(3,5);

  A(1,9) = C(1,3);
  A(2,9) = C(4,3);
  A(3,9) = C(6,3) + Sig(3,1);
  A(4,9) = C(4,3);
  A(5,9) = C(2,3);
  A(6,9) = C(5,3) + Sig(3,2);
  A(7,9) = C(6,3);
  A(8,9) = C(5,3);
  A(9,9) = C(3,3) + Sig(3,3);
}


bool NonlinearElasticityFbar::evalInt (LocalIntegral& elmInt,
				       const FiniteElement& fe,
				       const TimeDomain& prm,
				       const Vec3& X) const
{
#if INT_DEBUG > 1
  std::cout <<"NonlinearElasticityFbar::u(int) = "<< fe.u;
  if (nsd > 1) std::cout <<" "<< fe.v;
  if (nsd > 2) std::cout <<" "<< fe.w;
  std::cout <<"\nNonlinearElasticityFbar::dNdX ="<< fe.dNdX;
#endif

  FbarMats& fbar = static_cast<FbarMats&>(elmInt);

  // Evaluate the deformation gradient, F, at current configuration
  Matrix B, dNdx;
  Tensor F(nDF);
  SymmTensor E(nsd,axiSymmetry);
  if (!this->kinematics(fbar.vec.front(),fe.N,fe.dNdX,X.x,B,F,E))
    return false;

  double J, Jbar = 0.0;
  bool lHaveStrains = !E.isZero(1.0e-16);
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    if (nDF == nsd)
      Fi.fill(F.ptr());
    else
      for (unsigned short int i = 1; i <= nsd; i++)
	for (unsigned short int j = 1; j <= nsd; j++)
	  Fi(i,j) = F(i,j);

    J = Fi.inverse();
    if (axiSymmetry) J *= F(3,3);
    if (J == 0.0) return false;

    if (eKm || iS)
    {
      // Push-forward the basis function gradients to current configuration
      dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
#if INT_DEBUG > 0
      std::cout <<"NonlinearElasticityFbar::J = "<< J
		<<"\nNonlinearElasticityFbar::dNdx ="<< dNdx;
#endif
    }
  }
  else
  {
    // Initial state, no deformation yet
    F = J = 1.0;
    dNdx = fe.dNdX;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = (axiSymmetry ? 2.0*M_PI*X.x : 1.0)*fe.detJxW*J;
  double r = axiSymmetry ? X.x + fbar.vec.front().dot(fe.N,0,nsd) : 0.0;

  Vector M;
  Matrix dMdx;
  if (fbar.myVolData.size() == 1)
  {
    // Only one volumetric sampling point (linear element)
    Jbar = fbar.myVolData.front().J;
    dMdx = fbar.myVolData.front().dNdx;
    if (axiSymmetry) M = fbar.myVolData.front().Nr;
  }
  else if (!fbar.myVolData.empty())
  {
    // Evaluate the Lagrange polynomials extrapolating the volumetric
    // sampling points (assuming a regular distribution over the element)
    Vector Nbar;
    int pbar3 = nsd > 2 ? fbar.pbar : 0;
    if (!Lagrange::computeBasis(Nbar,
				fbar.pbar,fe.xi*fbar.scale,
				fbar.pbar,fe.eta*fbar.scale,
				pbar3,fe.zeta*fbar.scale)) return false;
#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityFbar::Nbar ="<< Nbar;
#endif

    // Compute modified deformation gradient determinant and basis function
    // gradients wrt. current configuration (spatial coordinates),
    // by extrapolating the volume sampling points
    dMdx.resize(dNdx.rows(),dNdx.cols(),true);
    if (axiSymmetry) M.resize(dNdx.rows(),true);
    for (size_t a = 0; a < Nbar.size(); a++)
    {
      Jbar += fbar.myVolData[a].J*Nbar[a];
      dMdx.add(fbar.myVolData[a].dNdx,fbar.myVolData[a].J*Nbar[a]);
      if (axiSymmetry)
	M.add(fbar.myVolData[a].Nr,fbar.myVolData[a].J*Nbar[a]);
    }
    dMdx.multiply(1.0/Jbar);
    if (axiSymmetry) M /= Jbar;
  }
  else
  {
    // No F-bar terms
    Jbar = J;
    dMdx = dNdx;
    if (axiSymmetry && r > 0.0)
      M = fe.N * (1.0/r);
  }

  // Compute modified deformation gradient, Fbar
  F *= pow(fabs(Jbar/J),1.0/(double)nDF);
#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityFbar::Jbar = "<< Jbar
	    <<"\nNonlinearElasticityFbar::dMdx ="<< dMdx
	    <<"NonlinearElasticityFbar::Fbar =\n"<< F;
#endif

  // Evaluate the constitutive relation (Jbar is dummy here)
  Matrix Cmat;
  SymmTensor Sig(nsd,axiSymmetry);
  if (!material->evaluate(Cmat,Sig,Jbar,fe.iGP,X,F,E,lHaveStrains,&prm))
    return false;

  // Multiply tangent moduli and stresses by integration point volume
  Cmat *= detJW;
  Sig  *= detJW;

  if (iS && lHaveStrains)
  {
    // Compute and accumulate contribution to internal force vector
    Vector& ES = fbar.b[iS-1];
    size_t a, d;
    unsigned short int i, j;
    for (a = d = 1; a <= dNdx.rows(); a++)
    {
      if (axiSymmetry && r > 0.0)
	ES(d) -= fe.N(a)*Sig(3,3)/r;
      for (i = 1; i <= nsd; i++, d++)
      {
	double Sint = 0.0;
	for (j = 1; j <= nsd; j++)
	  Sint += dNdx(a,j)*Sig(i,j);
	ES(d) -= Sint;
      }
    }
  }

  if (eKm)
  {
    // Compute standard and modified discrete spatial gradient operators
    Matrix G, Gbar;
    getGradOperator(axiSymmetry ? r : -1.0, fe.N, dNdx, G);
    getGradOperator(axiSymmetry ? 1.0 : -1.0, M, dMdx, Gbar);

    // Convert the spatial constitutive tensor to first elasticity tensor, A
    Matrix A;
    if (axiSymmetry)
      getAmatAx(Cmat,Sig,A);
    else if (nsd == 2)
      getAmat2D(Cmat,Sig,A);
    else
      getAmat3D(Cmat,Sig,A);

    // Compute the fourth-order tensor Q
    Matrix Q(A.rows(),A.cols());
    unsigned short int i, j, k, l;
    for (i = k = 1; k <= A.rows(); i++)
      for (j = 1; j <= nsd && k <= A.rows(); j++, k++)
      {
	Q(k,1) = (double)(1-nDF)*(i > nsd ? Sig(3,3) : Sig(i,j));
	for (l = 1; l <= nsd*nsd; l += nsd+1)
	  Q(k,1) += A(k,l);
	if (axiSymmetry)
	  Q(k,1) += A(k,5);
	Q(k,1) /= (double)nDF;
	for (l = nsd+2; l <= nsd*nsd; l += nsd+1)
	  Q(k,l) = Q(k,1);
	if (axiSymmetry)
	  Q(k,5) = Q(k,1);
      }

#if INT_DEBUG > 0
    std::cout <<"NonlinearElasticityFbar::G ="<< G
	      <<"NonlinearElasticityFbar::Gbar ="<< Gbar
	      <<"NonlinearElasticityFbar::A ="<< A
	      <<"NonlinearElasticityFbar::Q ="<< Q;
#endif

    // Compute and accumulate contribution to tangent stiffness matrix.

    // First, standard (material and geometric) tangent stiffness terms
    Matrix CB; CB.multiply(A,G);
    fbar.A[eKm-1].multiply(G,CB,true,false,true); // EK += G^T * A*G

    // Then, modify the tangent stiffness for the F-bar terms
    Gbar -= G; CB.multiply(Q,Gbar);
    fbar.A[eKm-1].multiply(G,CB,true,false,true); // EK += G^T * Q*Gbar
  }

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(fbar.A[eM-1],fe.N,X,detJW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(fbar.b[eS-1],fe.N,X,detJW);

  return true;
}


NormBase* NonlinearElasticityFbar::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormFbar(*const_cast<NonlinearElasticityFbar*>(this));
}


LocalIntegral* ElasticityNormFbar::getLocalIntegral (size_t nen, size_t iEl,
						     bool neumann) const
{
  LocalIntegral* result = this->NormBase::getLocalIntegral(nen,iEl,neumann);
  if (result) return new FbarNorm(static_cast<ElmNorm&>(*result));

  // Element norms are not requested, so allocate an internal object instead
  // that will delete itself when invoking its destruct method.
  return new FbarNorm(this->getNoFields());
}


bool ElasticityNormFbar::initElement (const std::vector<int>& MNPC,
				      const Vec3&, size_t nPt,
				      LocalIntegral& elmInt)
{
  NonlinearElasticityFbar& p = static_cast<NonlinearElasticityFbar&>(myProblem);
  FbarNorm& fbar = static_cast<FbarNorm&>(elmInt);

  return fbar.init(p.nsd,nPt) && this->NormBase::initElement(MNPC,*fbar.myNorm);
}


bool ElasticityNormFbar::initElementBou (const std::vector<int>& MNPC,
					 LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,*static_cast<FbarNorm&>(elmInt).myNorm);
}


bool ElasticityNormFbar::reducedInt (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X) const
{
  return myProblem.reducedInt(elmInt,fe,X);
}


bool ElasticityNormFbar::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const TimeDomain& prm,
				  const Vec3& X) const
{
  size_t nsd = fe.dNdX.cols();
  NonlinearElasticityFbar& p = static_cast<NonlinearElasticityFbar&>(myProblem);
  FbarNorm& fbar = static_cast<FbarNorm&>(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Matrix B;
  Tensor F(p.nDF);
  SymmTensor E(p.nDF);
  if (!p.kinematics(fbar.myNorm->vec.front(),fe.N,fe.dNdX,X.x,B,F,E))
    return false;

  double Jbar = 0.0;
  if (fbar.myVolData.size() == 1)
    Jbar = fbar.myVolData.front().J;
  else if (!fbar.myVolData.empty())
  {
    // Evaluate the Lagrange polynomials extrapolating the volumetric
    // sampling points (assuming a regular distribution over the element)
    RealArray Nbar;
    int pbar3 = nsd > 2 ? fbar.pbar : 0;
    if (!Lagrange::computeBasis(Nbar,
				fbar.pbar,fe.xi*fbar.scale,
				fbar.pbar,fe.eta*fbar.scale,
				pbar3,fe.zeta*fbar.scale)) return false;

    // Compute modified deformation gradient determinant
    // by extrapolating the volume sampling points
    for (size_t a = 0; a < Nbar.size(); a++)
      Jbar += fbar.myVolData[a].J*Nbar[a];
  }
  else
    Jbar = F.det();

  // Compute modified deformation gradient, Fbar
  Tensor Fbar(F);
  Fbar *= pow(fabs(Jbar/F.det()),1.0/(double)p.nDF);

  // Compute the strain energy density, U(E) = Int_E (S:Eps) dEps
  // and the Cauchy stress tensor, sigma
  Matrix Cmat; double U = 0.0;
  SymmTensor sigma(nsd, p.isAxiSymmetric() || p.material->isPlaneStrain());
  if (!p.material->evaluate(Cmat,sigma,U,fe.iGP,X,Fbar,E,3,&prm,&F))
    return false;

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = p.isAxiSymmetric() ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Integrate the norms
  return ElasticityNormUL::evalInt(*fbar.myNorm,sigma,U,F.det(),detJW);
}


bool ElasticityNormFbar::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& normal) const
{
  return this->ElasticityNormUL::evalBou(*static_cast<FbarNorm&>(elmInt).myNorm,
					 fe,X,normal);
}
