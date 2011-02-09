// $Id: StabilizedNavierStokes.C,v 1.5 2011-02-08 12:43:49 rho Exp $
//==============================================================================
//!
//! \file StabilizedNavierStokes.C
//!
//! \date Mar 17 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for stabilized Navier-Stokes problems.
//!
//==============================================================================

#include "StabilizedNavierStokes.h"
#include "Utilities.h"
#include "Vec3Oper.h"


StabilizedNavierStokes::StabilizedNavierStokes(short int n, ProblemFormulation form, int itg) 
  : StabilizedStokes(n,form,itg)
{
  mu  = 1.0e-3;
  rho = 1.0;
}


StabilizedNavierStokes::~StabilizedNavierStokes()
{
}


bool StabilizedNavierStokes::evalInt(LocalIntegral*& elmInt, double detJW,
			       const Vector& N, const Matrix& dNdX,
			       const Vec3& X) const
{
  int    i, j;
  int    k, l;
  double div, laplace, conv;
  Vector vel(3);

  if (eM) {
    Matrix& EM = *eM;
    Vector& EV = *eVs[0];
    
    // Velocity in integration point
    vel.fill(0.0);
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	vel(k) += EV((i-1)*nf + k)*N(i);

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += vel(k)*dNdX(j,k);
	conv *= rho*N(i)*detJW;

	for (k = 1;k <= nsd;k++)
	  EM((i-1)*nf+k,(j-1)*nf+k) += conv;
      }

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
    
    if (formulation == SIM::STRESS) 
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

  return getIntegralResult(elmInt);
}


bool StabilizedNavierStokes::evalInt (LocalIntegral*& elmInt, double detJW,
                                const Vector& N, const Matrix& dNdX,
                                const Matrix3D& d2NdX2, 
				const Vec3& X, double h) const
{
  int    i, j;
  int    k, l;
  double div, laplace, conv;
  Vector vel(3);

  const real delta = 0.001*h*h*detJW;

  if (eM) {
    Matrix& EM = *eM;
    Vector& EV = *eVs[0];

    // Velocity in integration point
    vel.fill(0.0);
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	vel(k) += EV((i-1)*nf + k)*N(i);

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += vel(k)*dNdX(j,k);
	conv *= rho*N(i)*detJW;

	for (k = 1;k <= nsd;k++)
	  EM((i-1)*nf+k,(j-1)*nf+k) += conv;
      }
    

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
    
    if (formulation == SIM::STRESS) 
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
    
    if (formulation == SIM::STRESS)
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

  return getIntegralResult(elmInt);
}


bool StabilizedNavierStokes::evalInt(LocalIntegral*& elmInt, double detJW,
			       const Vector& N, const Matrix& dNdX,
			       const Vector& Navg, const Vec3& X) const
{
  int    i, j;
  int    k, l;
  double div, laplace, conv;
  Vector vel(3);

  if (eM) {
    Matrix& EM = *eM;
    Vector& EV = *eVs[0];
    
    // Velocity in integration point
    vel.fill(0.0);
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	vel(k) += EV((i-1)*nf + k)*N(i);

    //Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += vel(k)*dNdX(j,k);
	conv *= rho*detJW;
	
	for (k = 1;k <= nsd;k++)
	  EM((i-1)*nf+k,(j-1)*nf+k) += conv*N(i);
      }
    
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
    
    if (formulation == SIM::STRESS) 
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

  return getIntegralResult(elmInt);
}


bool StabilizedNavierStokes::evalBou (LocalIntegral*& elmInt, double detJW,
                                const Vector& N, const Matrix& dNdX,
                                const Vec3& X, const Vec3& normal) const
{
  if (!eS || !tracFld)
  {
    std::cerr <<" *** StabilizedNavierStokes::evalBou: Zero pointers."<< std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec3 T = (*tracFld)(X,normal);

  // Store the traction value for vizualization
  if (!T.isZero())
    tracVal.insert(std::make_pair(X,T));

  // Integrate the force vector
  Vector& ES = *eS;
  for (size_t a = 1; a <= N.size(); a++)
    for (short int d = 1; d <= nsd; d++)
      ES(nf*(a-1)+d) += T[d-1]*N(a)*detJW;

  return getIntegralResult(elmInt);
}
