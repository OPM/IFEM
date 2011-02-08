//==============================================================================
//!
//! \file ChorinVelPredBE.C
//!
//! \date Sep 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for velocity prediction in Chorin's method
//!
//==============================================================================

#include "ChorinVelPredBE.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Utilities.h"


ChorinVelPredBE::ChorinVelPredBE(short int n, ProblemFormulation form,
				 int itg, bool incPress, bool mixed)
  : ChorinVelPred(n,form,itg,incPress,mixed)
{
  int i;

  // Resize velocity element solution vectors
  int nUsols = 2;
  eVs.resize(nUsols);
  for (i = 0;i < nUsols;i++)
    eVs[i] = new Vector();
  primsol.resize(nUsols);

  // Resize pressure element solution vectors
  int nPsols = 2;
  ePs.resize(nPsols);
  for (i = 0;i < nPsols;i++)
    ePs[i] = new Vector();
  psol.resize(nPsols);
}


ChorinVelPredBE::~ChorinVelPredBE() {}


bool ChorinVelPredBE::evalInt(LocalIntegral*& elmInt,
			    const TimeDomain& time, double detJW,
			    const Vector& N, const Matrix& dNdX,
			    const Vec3& X) const
{
  int i, j;
  int k, l;
  double dtInv, temp, mass, P;
  double conv, diff, laplace;
  Vector U(nsd), UP(nsd), dPdX(nsd);

  // Inverse of timestep size
  if (time.dt > 0.0)
    dtInv = 1.0/time.dt;
  else
    dtInv = 0.0;

  // Element velocity vectors
  Vector& EV  = *eVs[0];
  Vector& EVP = *eVs[0]; 

  // Velocity in integration point
  U.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      U(k) += EV((i-1)*nsd + k)*N(i);

  // Velocity in integration point at previous timestep
  UP.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      UP(k) += EVP((i-1)*nsd + k)*N(i);

  // Element pressure vector
  Vector& EP = *ePs[0];

  // Pressure in integration point
  P = 0.0;
  for (i = 1;i <= N.size();i++)
    P += EP(i)*N(i);
  
  // Pressure gradient in integration point
//   dPdX.fill(0.0);
//   for (i = 1;i <= N.size();i++)
//     for (k = 1;k <= nsd;k++)
//       dPdX(k) += EP(i)*dNdX(i,k);

  // Loop over basis functions

  // Lefthand size terms
  if (eM) {
    Matrix& EM = *eM;
    
    // Time derivative terms
    temp = rho*dtInv*detJW;
    for (i = 1;i <= N.size();i++) 
      for (j = 1;j <= N.size();j++) {
	mass = temp*N(i)*N(j);

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += mass;
      }

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += U(k)*dNdX(i,k);
	conv *= rho*N(j)*detJW;

	// Same convection contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += conv;
      }

    // Viscous terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= mu*detJW;

	// Same viscous contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += laplace;

	// Extra terms for stress-strain formulation
	if (formulation == STRESS)
	  for (k = 1;k <= nsd;k++)
	    for (l = 1;l <= nsd;l++)
	      EM((j-1)*nsd + l,(i-1)*nsd + k) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
      }
  }

  // Righthand side terms
  if (eS) {
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    // Vector fb(nsd);
//     this->getBodyForce(X,fb);
//     fb *= detJW;

    for (i = 1;i <= N.size();i++) {
      // Time derivative terms
      mass = rho*dtInv*N(i)*detJW;
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + k) += mass*UP(k); 

      // Body forces
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + k) += fb(k)*N(i);
    }

    // Only for incremental pressure
    if (incPressure) 
      for (i = 1;i <= N.size();i++)
	for (k = 1;k <= nsd;k++)
	  ES((i-1)*nsd + k) += P*dNdX(i,k)*detJW;
    //ES((i-1)*nsd + k) -= dPdX(k)*N(i)*detJW;
  }

  return getIntegralResult(elmInt);
}


bool ChorinVelPredBE::evalInt(LocalIntegral*& elmInt,
			    const TimeDomain& time, double detJW,
			    const Vector& N, const Matrix& dNdX,
			    const Matrix3D& d2NdX2, const Vec3& X,
			    double h) const
{
  int i, j;
  int k, l, m;
  double P, Uabs, delta;
  double dtInv, dUdt, temp;
  double conv, convI, convJ, diff, laplace;
  Vector U(nsd), UP(nsd), dPdX(nsd), res(nsd);
  Matrix   dUdX(nsd,nsd);
  Matrix3D d2UdX2(nsd,nsd,nsd);

  // Inverse of timestep size
  if (time.dt > 0.0)
    dtInv = 1.0/time.dt;
  else
    dtInv = 0.0;

  // Element velocity vectors
  Vector& EV  = *eVs[0];
  Vector& EVP = *eVs[0];
  
  // Velocity in integration point
  U.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      U(k) += EV((i-1)*nsd + k)*N(i);
  
  // Velocity magnitude
  Uabs = U.norm2();
  
  // Previous velocity in integration point
  UP.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      UP(k) += EVP((i-1)*nsd + k)*N(i);
  
  // Element pressure vector
  Vector& EP = *ePs[0];
  
  // Pressure in integration point
  P = 0.0;
  for (i = 1;i <= N.size();i++)
    P += EP(i)*N(i);

  // Velocity gradient in integration point
  dUdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	dUdX(k,l) += dNdX(i,l)*EV((i-1)*nsd + k);

  // Pressure gradient in integration point
  dPdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdX(k) += dNdX(i,k)*EP(i);

  // Hessian for velocity
  d2UdX2.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	for (m = 1;m <= nsd;m++)
	  d2UdX2(k,m,l) += d2NdX2(i,m,l)*EV((i-1)*nsd + k);

  // Stabilization parameters
  if (mu < rho*Uabs*h) {
    if (time.dt == 0)
      delta = 0.5/sqrt(1.0/(h*h) + (Uabs*Uabs)/(h*h));
    else
      delta = 0.5/sqrt(dtInv*dtInv + (Uabs*Uabs)/(h*h));
  }
  delta *= detJW;

  // Loop over basis functions

  // Lefthand size terms
  if (eM) {
    Matrix& EM = *eM;
    
    // Time derivative terms
    temp = rho*dtInv*detJW;
    for (i = 1;i <= N.size();i++) 
      for (j = 1;j <= N.size();j++) {
	dUdt = temp*N(i)*N(j);

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += dUdt;
      }

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += U(k)*dNdX(i,k);
	conv *= rho*N(j)*detJW;

	// Same convection contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += conv;
      }

    // Viscous terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= mu*detJW;
	
	// Same viscous contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += laplace;

	// Extra terms for stress-strain formulation
	if (formulation == STRESS)
	  for (k = 1;k <= nsd;k++)
	    for (l = 1;l <= nsd;l++)
	      EM((j-1)*nsd + l,(i-1)*nsd + k) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
      }

    // SUPG stabilization terms
    for (i = 1;i <= N.size();i++) {
      // Convection for solution basis functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
	convI += U(k)*dNdX(i,k);
      convI *= rho;

      for (j = 1;j <= N.size();j++) {
	// Convection for test functions
	convJ = 0.0;
	for (k = 1;k <= nsd;k++)
	  convJ += U(k)*dNdX(j,k);
	convJ *= rho*delta;

	// Convection terms
	for (k = 1;k <= nsd;k++) 
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += convI*convJ;

	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= mu;

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) -= diff*convJ;
      }
    }
  }

  // Righthand side terms
  if (eS) {
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X);

    for (i = 1;i <= N.size();i++) {
      // Time derivative terms
      dUdt = rho*N(i)*dtInv*detJW;
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + k) += dUdt*UP(k); 

      // Body forces
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + l) += fb(k)*N(i);

      // Convection for test functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
      	convI += U(k)*dNdX(i,k);
      convI *= rho*delta;

      // Body force
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + l) += fb(k)*(N(i)*detJW + convI);

      // Pressure term
      if (incPressure)
	for (k = 1;k <= nsd;k++)
	  ES((i-1)*nsd + k) += P*dNdX(i,k)*detJW - dPdX(k)*convI;
	  //	  ES((i-1)*nsd + k) -= dPdX(k)*(N(i)*detJW + convI);
    }      
  }
  
  return getIntegralResult(elmInt);
}


bool ChorinVelPredBE::evalInt(LocalIntegral*& elmInt, 
			      const TimeDomain& time, double detJW,
			      const Vector& N1, const Vector& N2,
			      const Matrix& dN1dX, const Matrix& dN2dX,
			      const Vec3& X) const
{
  int i, j;
  int k, l;
  double dtInv, temp, mass, P;
  double conv, diff, laplace;
  Vector U(nsd), UP(nsd), dPdX(nsd);

  // Inverse of timestep size
  if (time.dt > 0.0)
    dtInv = 1.0/time.dt;
  else
    dtInv = 0.0;

  // Element velocity vectors
  Vector& EV  = *eVs[0];
  Vector& EVP = *eVs[0]; 

  // Velocity in integration point
  U.fill(0.0);
  for (i = 1;i <= N1.size();i++)
    for (k = 1;k <= nsd;k++)
      U(k) += EV((i-1)*nsd + k)*N1(i);

  // Velocity in integration point at previous timestep
  UP.fill(0.0);
  for (i = 1;i <= N1.size();i++)
    for (k = 1;k <= nsd;k++)
      UP(k) += EVP((i-1)*nsd + k)*N1(i);

  // Element pressure vector
  Vector& EP = *ePs[0];

  // Pressure in integration point
  P = 0.0;
  for (i = 1;i <= N2.size();i++)
    P += EP(i)*N2(i);

  // Loop over basis functions

  // Lefthand size terms
  if (eM) {
    Matrix& EM = *eM;
    
    // Time derivative terms
    temp = rho*dtInv*detJW;
    for (i = 1;i <= N1.size();i++) 
      for (j = 1;j <= N1.size();j++) {
	mass = temp*N1(i)*N1(j);

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += mass;
      }

    // Convection terms
    for (i = 1;i <= N1.size();i++)
      for (j = 1;j <= N1.size();j++) {
	// Sum convection for each direction
	conv = 0.0;
	for (k = 1;k <= nsd;k++)
	  conv += U(k)*dN1dX(i,k);
	conv *= rho*N1(j)*detJW;

	// Same convection contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += conv;
      }

    // Viscous terms
    for (i = 1;i <= N1.size();i++)
      for (j = 1;j <= N1.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dN1dX(i,k)*dN1dX(j,k);
	laplace *= mu*detJW;

	// Same viscous contribution to each velocity component
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd + k,(i-1)*nsd + k) += laplace;

	// Extra terms for stress-strain formulation
	if (formulation == STRESS)
	  for (k = 1;k <= nsd;k++)
	    for (l = 1;l <= nsd;l++)
	      EM((j-1)*nsd + l,(i-1)*nsd + k) += mu*dN1dX(i,k)*dN1dX(j,l)*detJW;
      }
  }

  // Righthand side terms
  if (eS) {
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;

    for (i = 1;i <= N1.size();i++) {
      // Time derivative terms
      mass = rho*dtInv*N1(i)*detJW;
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + k) += mass*UP(k); 

      // Body forces
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nsd + k) += fb(k)*N1(i);
    }

    // Only for incremental pressure
    if (incPressure) 
      for (i = 1;i <= N1.size();i++)
	for (k = 1;k <= nsd;k++)
	  ES((i-1)*nsd + k) += P*dN1dX(i,k)*detJW;
  }

  return getIntegralResult(elmInt); 
}


bool ChorinVelPredBE::evalBou(LocalIntegral*& elmInt, 
			      const TimeDomain& time, double detJW,
			      const Vector& N, const Matrix& dNdX,
			      const Vec3& X, const Vec3& normal) const
{
  if (incPressure) {
    if (!tracFld) {
      std::cerr <<" *** ChorinVelPredBE::evalBou: No tractions."<< std::endl;
      return false;
    }
    else if (!eS) {
      std::cerr <<" *** ChorinVelPredBE::evalBou: No load vector."<< std::endl;
      return false;
    }
    
    // Evaluate the surface traction
    Vec3 T = (*tracFld)(X,normal);
    
    // Store the traction value for vizualization
    if (!T.isZero()) tracVal[X] = T;
    
    // Integrate the force vector
    Vector& ES = *eS;
    for (size_t a = 1; a <= N.size(); a++)
      for (short int d = 1; d <= nsd; d++)
	ES(nsd*(a-1)+d) -= T[d-1]*N(a)*detJW;
  }
  
  return getIntegralResult(elmInt);
}

