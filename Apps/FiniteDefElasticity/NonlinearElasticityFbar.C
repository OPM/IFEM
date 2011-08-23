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
#include "Lagrange.h"
#include "ElmNorm.h"
#include "TimeDomain.h"
#include "Tensor.h"
#include "GaussQuadrature.h"


NonlinearElasticityFbar::NonlinearElasticityFbar (unsigned short int n, int nvp)
  : NonlinearElasticityUL(n), npt1(nvp)
{
  iP = 0;
  pbar = 0;
  scale = 1.0;
}


void NonlinearElasticityFbar::print (std::ostream& os) const
{
  std::cout <<"NonlinearElasticityFbar: F-bar formulation, "<< npt1
	    <<" volumetric points in each direction."<< std::endl;

  this->NonlinearElasticityUL::print(os);
}


bool NonlinearElasticityFbar::initElement (const std::vector<int>& MNPC,
                                           const Vec3&, size_t nPt)
{
  iP = 0;
  pbar = (int)pow((double)nPt,1.0/(double)nsd);
  if (pow((double)pbar,(double)nsd) == (double)nPt)
    myVolData.resize(nPt);
  else
  {
    pbar = 0;
    std::cerr <<" *** NonlinearElasticityFbar::initElement: Invalid element, "
	      << nPt <<" volumetric sampling points specified."<< std::endl;
    return false;
  }
  scale = pbar > 1 ? 1.0/GaussQuadrature::getCoord(pbar)[pbar-1] : 1.0;
#if INT_DEBUG > 0
  std::cout <<"\n\n *** Entering NonlinearElasticityFbar::initElement";
  std::cout <<", nPt = "<< nPt <<", pbar = "<< pbar
	    <<", scale = "<< scale << std::endl;
#endif

  return this->NonlinearElasticityUL::initElement(MNPC);
}


bool NonlinearElasticityFbar::reducedInt (const FiniteElement& fe) const
{
#if INT_DEBUG > 1
  std::cout <<"NonlinearElasticityFbar::u(red) = "<< fe.u;
  if (nsd > 1) std::cout <<" "<< fe.v;
  if (nsd > 2) std::cout <<" "<< fe.w;
  std::cout <<" (iP="<< iP+1 <<")\n"
	    <<"NonlinearElasticityFbar::dNdX ="<< fe.dNdX;
#endif

  VolPtData& ptData = myVolData[iP++];

  // Evaluate the deformation gradient, F, at current configuration
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(fe.dNdX,F,E))
    return false;

  if (E.isZero(1.0e-16))
  {
    // Initial state, no deformation yet
    ptData.J = 1.0;
    ptData.dNdx = fe.dNdX;
  }
  else
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    Fi.fill(F.ptr());
    ptData.J = Fi.inverse();
    if (ptData.J == 0.0)
      return false;

    // Push-forward the basis function gradients to current configuration
    ptData.dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
#ifdef INT_DEBUG
    std::cout <<"NonlinearElasticityFbar::J = "<< ptData.J
	      <<"\nNonlinearElasticityFbar::dNdx ="<< ptData.dNdx;
#endif
  }
#ifdef INT_DEBUG
  if (iP == myVolData.size())
    std::cout <<"NonlinearElasticityFbar: Volumetric sampling points completed."
	      <<"\n"<< std::endl;
#endif

  return true;
}


/*!
  \brief Computes the discrete gradient operator, \b G, for continuum elements.
*/

static void getGradOperator (const Matrix& dNdx, Matrix& G)
{
  size_t nen = dNdx.rows();
  size_t nsd = dNdx.cols();
  G.resize(nsd*nsd,nsd*nen);

  size_t a, i, j, k, l = 0;
  for (a = 1; a <= nen; a++, l += nsd)
    for (i = k = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++, k++)
	G(k,l+j) = dNdx(a,i);
}


/*!
  \brief Computes the spatial tangent modulus matrix, \b A, for 2D elements.
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
  \brief Computes the spatial tangent modulus matrix, \b A, for 3D elements.
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


bool NonlinearElasticityFbar::evalInt (LocalIntegral*& elmInt,
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

  // Evaluate the deformation gradient, F, at current configuration
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!this->kinematics(fe.dNdX,F,E))
    return false;

  double J, Jbar = 0.0;
  bool lHaveStrains = !E.isZero(1.0e-16);
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    Fi.fill(F.ptr());
    J = Fi.inverse();
    if (J == 0.0) return false;

    if (eKm || iS)
    {
      // Push-forward the basis function gradients to current configuration
      dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
#ifdef INT_DEBUG
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

  if (myVolData.size() == 1)
  {
    // Only one volumetric sampling point (linear element)
    Jbar = myVolData.front().J;
    dMdx = myVolData.front().dNdx;
  }
  else if (!myVolData.empty())
  {
    // Evaluate the Lagrange polynomials extrapolating the volumetric
    // sampling points (assuming a regular distribution over the element)
    RealArray Nbar;
    int pbar3 = nsd > 2 ? pbar : 0;
    if (Lagrange::computeBasis(Nbar,pbar,fe.xi*scale,pbar,fe.eta*scale,
			       pbar3,fe.zeta*scale)) return false;

#ifdef INT_DEBUG
    std::cout <<"NonlinearElasticityFbar::Nbar =";
    for (size_t i = 0; i < Nbar.size(); i++) std::cout <<" "<< Nbar[i];
    std::cout << std::endl;
#endif

    // Compute modified deformation gradient determinant and basis function
    // gradients wrt. current configuration (spatial coordinates),
    // by extrapolating the volume sampling points
    dMdx.resize(dNdx.rows(),dNdx.cols(),true);
    for (size_t a = 0; a < Nbar.size(); a++)
    {
      Jbar += myVolData[a].J*Nbar[a];
      dMdx.add(myVolData[a].dNdx,myVolData[a].J*Nbar[a]);
    }
    dMdx.multiply(1.0/Jbar);
  }
  else
  {
    // No F-bar terms
    Jbar = J;
    dMdx = dNdx;
  }

  // Compute modified deformation gradient, Fbar
  F *= pow(fabs(Jbar/J),1.0/(double)nsd);
#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityFbar::Jbar = "<< Jbar
	    <<"\nNonlinearElasticityFbar::dMdx ="<< dMdx
	    <<"NonlinearElasticityFbar::Fbar =\n"<< F;
#endif

  // Evaluate the constitutive relation (Jbar is dummy here)
  SymmTensor Sig(nsd);
  if (!material->evaluate(Cmat,Sig,Jbar,X,F,E,lHaveStrains,&prm))
    return false;

  // Multiply tangent moduli and stresses by integration point volume
  Cmat *= fe.detJxW*J;
  Sig  *= fe.detJxW*J;

  if (iS && lHaveStrains)
  {
    // Compute and accumulate contribution to internal force vector
    Vector& ES = *iS;
    size_t a, d;
    unsigned short int i, j;
    for (a = d = 1; a <= dNdx.rows(); a++)
      for (i = 1; i <= nsd; i++, d++)
      {
	double Sint = 0.0;
	for (j = 1; j <= nsd; j++)
	  Sint += dNdx(a,j)*Sig(i,j);
	ES(d) -= Sint;
      }
  }

  if (eKm)
  {
    // Compute standard and mixed discrete spatial gradient operators
    Matrix G, Gbar;
    getGradOperator(dNdx,G);
    getGradOperator(dMdx,Gbar);

    // Convert the spatial constitutive tensor to first elasticity tensor, A
    Matrix A(nsd*nsd,nsd*nsd), Q(nsd*nsd,nsd*nsd);
    if (nsd == 2)
      getAmat2D(Cmat,Sig,A);
    else
      getAmat3D(Cmat,Sig,A);

    // Compute the fourth-order tensor Q
    unsigned short int i, j, k, l;
    for (i = k = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++, k++)
      {
	Q(k,1) = (double)(1-nsd)*Sig(i,j);
	for (l = 1; l <= A.cols(); l += nsd+1)
	  Q(k,1) += A(k,l);
	Q(k,1) /= (double)nsd;
	for (l = nsd+2; l <= Q.cols(); l += nsd+1)
	  Q(k,l) += Q(k,1);
      }

#ifdef INT_DEBUG
  std::cout <<"NonlinearElasticityFbar::G ="<< G
	    <<"NonlinearElasticityFbar::Gbar ="<< Gbar
	    <<"NonlinearElasticityFbar::A ="<< A
	    <<"NonlinearElasticityFbar::Q ="<< Q;
#endif

    // Compute and accumulate contribution to tangent stiffness matrix
    // First, standard (material and geometric) tangent stiffness terms
    eKm->multiply(G,CB.multiply(A,G),true,false,true); // EK += G^T * A*G

    // Modify the tangent stiffness for the F-bar terms
    Gbar -= G;
    eKm->multiply(G,CB.multiply(Q,Gbar),true,false,true); // EK += G^T * Q*Gbar
  }

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(*eM,fe.N,X,J*fe.detJxW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(*eS,fe.N,X,J*fe.detJxW);

  return this->getIntegralResult(elmInt);
}



NormBase* NonlinearElasticityFbar::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormFbar(*const_cast<NonlinearElasticityFbar*>(this));
}


bool ElasticityNormFbar::reducedInt (const FiniteElement& fe) const
{
  return myProblem.reducedInt(fe);
}


bool ElasticityNormFbar::evalInt (LocalIntegral*& elmInt,
				  const FiniteElement& fe,
				  const TimeDomain& prm,
				  const Vec3& X) const
{
  NonlinearElasticityFbar& p = static_cast<NonlinearElasticityFbar&>(myProblem);

  ElmNorm& pnorm = NormBase::getElmNormBuffer(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  size_t nsd = fe.dNdX.cols();
  Tensor F(nsd);
  SymmTensor E(nsd);
  if (!p.kinematics(fe.dNdX,F,E))
    return false;

  double Jbar = 0.0;
  if (p.myVolData.size() == 1)
    Jbar = p.myVolData.front().J;
  else if (!p.myVolData.empty())
  {
    // Evaluate the Lagrange polynomials extrapolating the volumetric
    // sampling points (assuming a regular distribution over the element)
    RealArray Nbar;
    int pbar3 = nsd > 2 ? p.pbar : 0;
    if (!Lagrange::computeBasis(Nbar,p.pbar,fe.xi*p.scale,p.pbar,fe.eta*p.scale,
				pbar3,fe.zeta*p.scale)) return false;

    // Compute modified deformation gradient determinant
    // by extrapolating the volume sampling points
    for (size_t a = 0; a < Nbar.size(); a++)
      Jbar += p.myVolData[a].J*Nbar[a];
  }
  else
    Jbar = F.det();

  // Compute modified deformation gradient, Fbar
  Tensor Fbar(F);
  Fbar *= pow(fabs(Jbar/F.det()),1.0/(double)nsd);

  // Compute the strain energy density, U(E) = Int_E (Sig:Eps) dEps
  double U = 0.0;
  SymmTensor S(nsd,p.material->isPlaneStrain());
  if (!p.material->evaluate(p.Cmat,S,U,X,Fbar,E,3,&prm,&F))
    return false;

  // Integrate the energy norm a(u^h,u^h) = Int_Omega0 U(E) dV0
  pnorm[0] += U*fe.detJxW;
  // Integrate the L2-norm ||Sig|| = Int_Omega0 Sig:Sig dV0
  pnorm[2] += S.L2norm(false)*fe.detJxW;
  // Integrate the von Mises stress norm
  pnorm[3] += S.vonMises(false)*fe.detJxW;

  return true;
}
