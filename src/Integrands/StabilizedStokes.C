// $Id: StabilizedStokes.C,v 1.6 2010-12-30 17:38:19 kmo Exp $
//==============================================================================
//!
//! \file StabilizedStokes.C
//!
//! \date Feb 11 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for stabilized Stokes problems.
//!
//==============================================================================

#include "StabilizedStokes.h"
#include "Utilities.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "VTF.h"


StabilizedStokes::StabilizedStokes (short int n,
				    Stokes::ProblemFormulation form, int itg)
  : Stokes(n,(SIM::Formulation)form,itg)
{
  // Resize element solution vectors
  eVs.resize(2);
  for (int i = 0;i < 2;i++)
    eVs[i] = new Vector();

  primsol.resize(2);
}


bool StabilizedStokes::evalInt (LocalIntegral*& elmInt, double detJW,
				const Vector& N, const Matrix& dNdX,
				const Vec3& X) const
{
  size_t    i, j;
  short int k, l;
  double    div, laplace;

  const int nf = nsd + 1;

  elmInt = &myMats;
  if (eM) {
    Matrix& EM = *eM;

    // Viscous terms
    for (i = 1; i <= N.size(); i++)
      for (j = 1; j <= N.size(); j++) {
	laplace = 0.0;
	for (k = 1; k <= nsd; k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= mu*detJW;
	for (k = 1; k <= nsd; k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    if (formulation == STRESS) 
      for (i = 1; i <= N.size(); i++)
 	for (j = 1; j <= N.size(); j++)
 	  for (k = 1; k <= nsd; k++)
 	    for (l = 1; l <= nsd; l++)
 	      EM(nf*(j-1)+k,nf*(i-1)+l) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
    
    // Pressure and divergence terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)
	for (k = 1;k <= nsd;k++) {
	  div = N(j)*dNdX(i,k)*detJW;
	  EM(nf*(i-1)+k,nf*j) -= div;
	  EM(nf*j,nf*(i-1)+k) -= div;
	}

    // Pressure stabilization
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) 
	EM(nf*j,nf*i) -= N(i)*N(j)*detJW;
  }
  
  if (eS) {
    // Integrate body force vector
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) += fb(k)*N(i);
  }	

  myMats.rhsOnly = false;
  return true;
}


bool StabilizedStokes::evalInt (LocalIntegral*& elmInt, double detJW,
                                const Vector& N, const Matrix& dNdX,
                                const Matrix3D& d2NdX2, 
				const Vec3& X, double h) const
{
  size_t    i, j;
  short int k, l;
  double    div, laplace;

  const int nf = nsd + 1;
  const real delta = 0.001*h*h*detJW;

  elmInt = &myMats;
  if (eM) {
    Matrix& EM = *eM;

    // Viscous terms
    for (i = 1; i <= N.size(); i++)
      for (j = 1; j <= N.size(); j++) {
	laplace = 0.0;
	for (k = 1; k <= nsd; k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= mu*detJW;
	for (k = 1; k <= nsd; k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    if (formulation == STRESS) 
      for (i = 1; i <= N.size(); i++)
 	for (j = 1; j <= N.size(); j++)
 	  for (k = 1; k <= nsd; k++)
 	    for (l = 1; l <= nsd; l++)
 	      EM(nf*(j-1)+k,nf*(i-1)+l) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
    
    // Pressure and divergence terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)
	for (k = 1;k <= nsd;k++) {
	  div = N(j)*dNdX(i,k)*detJW;
	  EM(nf*(i-1)+k,nf*j) -= div;
	  EM(nf*j,nf*(i-1)+k) -= div;
	}

    // Pressure stabilization
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	for (k = 1;k <= nsd;k++) {
  	  laplace = 0.0;
  	  for (l = 1;l <= nsd;l++)
  	    laplace += d2NdX2(i,l,l);
  	  laplace *= delta*mu;
	  
 	  EM(nf*j,nf*(i-1) + k) += laplace*dNdX(j,k);
	  EM(nf*j,nf*i) -= delta*dNdX(i,k)*dNdX(j,k);
	}
      }
    
    if (formulation == STRESS)
      for (i = 1;i <= N.size();i++)
 	for (j = 1;j <= N.size();j++)
 	  for (k = 1;k <= nsd;k++)
 	    for (l = 1;l <= nsd;l++)
 	      if (k != l) 
 		EM(nf*j,nf*(i-1)+l) += delta*mu*d2NdX2(i,k,l)*dNdX(j,k);
  }
  
  if (eS) {
    // Integrate body force vector
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) += fb(k)*N(i);

    // Pressure stabilization term
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
     	ES(nf*i) -= delta*fb(k)*dNdX(i,k);
  }	

  myMats.rhsOnly = false;
  return true;
}


bool StabilizedStokes::evalInt (LocalIntegral*& elmInt, double detJW,
				const Vector& N, const Matrix& dNdX,
				const Vector& Navg, const Vec3& X) const
{
  size_t    i, j;
  short int k, l;
  double    div, laplace;

  const int nf = nsd + 1;

  elmInt = &myMats;
  if (eM) {
    Matrix& EM = *eM;

    // Viscous terms
    for (i = 1; i <= N.size(); i++)
      for (j = 1; j <= N.size(); j++) {
	laplace = 0.0;
	for (k = 1; k <= nsd; k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= mu*detJW;
	for (k = 1; k <= nsd; k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    if (formulation == STRESS) 
      for (i = 1; i <= N.size(); i++)
 	for (j = 1; j <= N.size(); j++)
 	  for (k = 1; k <= nsd; k++)
 	    for (l = 1; l <= nsd; l++)
 	      EM(nf*(j-1)+k,nf*(i-1)+l) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
    
    // Pressure and divergence terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)
	for (k = 1;k <= nsd;k++) {
	  div = N(j)*dNdX(i,k)*detJW;
	  EM(nf*(i-1)+k,nf*j) -= div;
	  EM(nf*j,nf*(i-1)+k) -= div;
	}

    // Pressure stabilization
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)  
	EM(nf*j,nf*i) -= (N(i)-Navg(i))*(N(j)-Navg(j))*detJW;
  }
  
  if (eS) {
    // Integrate body force vector
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) += fb(k)*N(i);
  }	

  myMats.rhsOnly = false;
  return true;
}


