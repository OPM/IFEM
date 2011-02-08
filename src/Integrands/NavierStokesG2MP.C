// $Id: NavierStokesG2MP.C,v 1.2 2010-10-05 09:26:12 kmo Exp $
//==============================================================================
//!
//! \file NavierStokesG2MP.C
//!
//! \date August 11 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for G2 stabilized Navier-Stokes problems.
//! \details Time integration by the Crank-Nicholson method.
//!
//==============================================================================

#include "NavierStokesG2MP.h"
#include "TimeDomain.h"
#include "ElmNorm.h"


NavierStokesG2MP::NavierStokesG2MP(short int n, ProblemFormulation form, int itg)
  : NavierStokesG2(n,form,itg) {}


NavierStokesG2MP::~NavierStokesG2MP() {}


bool NavierStokesG2MP::evalInt(LocalIntegral*& elmInt, 
				     const TimeDomain& time, double detJW,
				     const Vector& N, const Matrix& dNdX,
				     const Matrix3D& d2NdX2, const Vec3& X,
				     double h) const
{
  int      i, j;
  int      k, l, m;
  real     temp1, temp2;
  real     diff, divU, divI, convI, convJ, laplace, U, P;
  real     Pb, Pn;
  real     delta1, delta2; 
  real     dtInv, dUdt;
  Vector   vel(3), bvel(3), conv(3), res(3), dPdX(3), dPdXn(3), dPdXb(3);
  Vector   pvel(3);
  Matrix   dUdX(3,3), dUdXn(3,3), dUdXb(3,3);
  Matrix3D d2UdX2(3,3,3), d2UdX2n(3,3,3), d2UdX2b(3,3,3);

  const int nf = nsd + 1;

  // Inverse of timestep size
  if (time.dt > 0.0)
    dtInv = 1.0/time.dt;
  else
    dtInv = 0.0;
  
  Vector& EV = *eVs[0];
  // Velocity in integration point
  vel.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      vel(k) += EV((i-1)*nf + k)*N(i);

  U = 0.0;
  for (k = 1;k <= nsd;k++)
    U += vel(k)*vel(k);
  U = sqrt(U);
  
  Vector& EVP = *eVs[1];
  // Previous velocity at integration point
  pvel.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++) 
      pvel(k) += EVP((i-1)*nf + k)*N(i);

  // Intermediate velocity at integration point
  bvel.fill(0.0);
  for (k = 1;k <= nsd;k++)
    bvel(k) = 0.5*(vel(k)+pvel(k));

  // Pressure in integration point
  P = 0.0;
  for (i = 1;i <= N.size();i++)
    P += EV(i*nf)*N(i);

  // Pressure in integration point at previous timestep
  Pn = 0.0;
  for (i = 1;i <= N.size();i++)
    Pn += EVP(i*nf)*N(i);

  // Pressure in integration point at intermediate timestep
  Pb = 0.5*(P+Pn);

  // Velocity gradient in integration point
  dUdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	dUdX(k,l) += dNdX(i,l)*EV((i-1)*nf+k);
  
  // Velocity gradient at previous timestep
  dUdXn.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	dUdXn(k,l) += dNdX(i,l)*EVP((i-1)*nf+k);

  // Velocity gradient at intermediate timestep
  for (k = 1;k <= nsd;k++)
    for (l = 1;l <= nsd;l++)
      dUdXb(k,l) = 0.5*(dUdX(k,l)+dUdXn(k,l));

  // Pressure gradient in integration point
  dPdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdX(k) += dNdX(i,k)*EV(i*nf);

  // Pressure gradient at previous timestep
  dPdXn.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdXn(k) += dNdX(i,k)*EVP(i*nf);

  // Pressure gradient at intermediate timestep
  for (k = 1;k <= nsd;k++)
    dPdXb(k) = 0.5*(dPdX(k) + dPdXn(k));

  // Hessian for velocity
  d2UdX2.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	for (m = 1;m <= nsd;m++)
	  d2UdX2(k,m,l) += d2NdX2(i,m,l)*EV((i-1)*nf + k);
  
  // Hessian for velocity at previous timestep
  d2UdX2n.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	for (m = 1;m <= nsd;m++)
	  d2UdX2n(k,m,l) += d2NdX2(i,m,l)*EVP((i-1)*nf + k);

  // Hessian for velocity at intermediate timestep
  for (k = 1;k <= nsd;k++)
    for (m = 1;m <= nsd;m++)
      for (l = 1;l <= nsd;l++)
	d2UdX2b(k,m,l) = 0.5*(d2UdX2(k,m,l) + d2UdX2n(k,m,l));

  // Stabilization parameters
  if (mu < rho*U*h) {
    // RUNAR
    if (time.dt == 0.0)
      delta1 = 0.5/sqrt(1.0/(h*h) + (U*U)/(h*h));
    else
      delta1 = 0.5/sqrt(dtInv*dtInv + (U*U)/(h*h));

    delta2 = h;
  }
  else 
    delta1 = delta2 = h*h;
  // RUNAR
  //delta1 *= 0.01*detJW;
  //delta2 *= 0.01*detJW;
  delta1 *= 0.001*detJW;
  delta2 *= 0.001*detJW;
  // RUNAR
  //delta1 = 0.0001*h*h*detJW;
  //delta2 = 0.0001*h*h*detJW;

  elmInt = &myMats;
  if (eM) {
    Matrix& EM = *eM;

    // Time derivative term lhs
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	dUdt = rho*N(i)*N(j)*dtInv*detJW;
	
	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nf + k, (i-1)*nf + k) += dUdt;
      }

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	convJ = 0.0;
	for (k = 1;k <= nsd;k++)
	  convJ += bvel(k)*dNdX(i,k);
	convJ *= 0.5*rho*N(j)*detJW;

	temp1 = 0.5*rho*N(i)*N(j)*detJW;

	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf+k,(i-1)*nf+k) += convJ;
	
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) += temp1*dUdXb(l,k);
	}
      }
     
    // Viscous terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= 0.5*mu*detJW;
	for (k = 1;k <= nsd;k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    temp1 = 0.5*mu*detJW;
    if (formulation == STRESS) 
      for (i = 1; i <= N.size(); i++)
 	for (j = 1; j <= N.size(); j++)
 	  for (k = 1; k <= nsd; k++)
 	    for (l = 1; l <= nsd; l++)
 	      EM(nf*(j-1)+l,nf*(i-1)+k) += temp1*dNdX(i,k)*dNdX(j,l);
    
    // Pressure and divergence terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)
	for (k = 1;k <= nsd;k++) {
	  divI = N(j)*dNdX(i,k)*detJW;
	  EM(nf*(i-1)+k,nf*j) -= 0.5*divI;
	  EM(nf*j,nf*(i-1)+k) -= 0.5*divI;
	}
    
    // Stabilization terms
    
    // Computation of non-linear residual
    res.fill(0.0);
    for (k = 1;k <= nsd;k++) {
      for (l = 1;l <= nsd;l++) 
	res(k)  += rho*bvel(l)*dUdXb(k,l)   - mu*d2UdX2b(k,l,l);

      res(k) += dPdXb(k);
      res(k) *= delta1;
    }

    // STABILIZATION TERMS
    for (i = 1;i <= N.size();i++) {
      // Convection for solution basis functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
	convI += bvel(k)*dNdX(i,k);
      convI *= rho;
      
      for (j = 1;j <= N.size();j++) {
	// Convection for test functions
	convJ = 0.0;
	for (k = 1;k <= nsd;k++)
	  convJ += bvel(k)*dNdX(j,k);
	convJ *= delta1*rho;
	
	//--- MOMENTUM STABILIZATION ---

	// Convection terms
	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf + k,(i-1)*nf + k) += 0.5*convI*convJ;
	  
	  for (l = 1;l <= nsd;l++) {
	    EM((j-1)*nf + l,(i-1)*nf + k) += 0.5*rho*res(l)*dNdX(j,k)*N(i);	 
	    EM((j-1)*nf + l,(i-1)*nf + k) += 0.5*convJ*rho*dUdXb(l,k)*N(i);	 
	  }
	}
	
	// Pressure terms
	for (l = 1;l <= nsd;l++) 
	  EM((j-1)*nf + l,i*nf) += 0.5*dNdX(i,l)*convJ;

	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= 0.5*mu;

	for (k = 1;k <= nsd;k++) 
	  EM((j-1)*nf + k,(i-1)*nf + k) -= diff*convJ;

 	// Continuity-continuity stabilization
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) -= 0.5*delta2*dNdX(i,k)*dNdX(j,l);

	//--- CONTINUITY STABILIZATION ---
	
	// Convection terms
	for (k = 1;k <= nsd;k++)
	  EM(j*nf,(i-1)*nf+k) += 0.5*delta1*convI*dNdX(j,k);
	
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM(j*nf,(i-1)*nf+k) += 0.5*delta1*rho*dUdXb(l,k)*N(i)*dNdX(j,l);
	
	// Pressure terms
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= 0.5*delta1;
	EM(j*nf,i*nf) += laplace;
	
	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= 0.5*delta1*mu;
	
	for (k = 1;k <= nsd;k++)
	  EM(j*nf,(i-1)*nf+k) -= diff*dNdX(j,k);
      }
    }
  }    

  if (eS) {
    // Integrate rhs vector
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    for (i = 1;i <= N.size();i++) {
      // --- Residual of Galerkin formulation ---

      // -- Momentum equations --

      // Time derivative terms
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf + k) += rho*(pvel(k)-vel(k))*N(i)*dtInv*detJW;

      // Body forces
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) += fb(k)*N(i);

      // Convection terms
      for (k = 1;k <= nsd;k++) {
	conv.fill(0.0);
	for (m = 1;m <= nsd;m++)
	  conv(k) += bvel(m)*dUdXb(k,m);

	ES((i-1)*nf + k) -= rho*conv(k)*N(i)*detJW;
      }

      // Pressure terms
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf + k) += Pb*dNdX(i,k)*detJW;

      // Viscous terms
      for (k = 1;k <= nsd;k++)
	for (m = 1;m <= nsd;m++)
	  ES((i-1)*nf + k) -= mu*dUdXb(k,m)*dNdX(i,m)*detJW;

      // -- Continuity equation --
      divU = 0.0;
      for (m = 1;m <= nsd;m++)
	divU += dUdXb(m,m);

      ES(i*nf) += divU*N(i)*detJW;

      //--- MOMENTUM STABILIZATION---

      // Convection for solution basis functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
	convI += vel(k)*dNdX(i,k);
      convI *= rho;

      // Convection terms
      for (k = 1;k <= nsd;k++) 
	ES((i-1)*nf+k) -= res(k)*convI; 

      // Continuity residual term
      for (k = 1;k <= nsd;k++)
      	ES((i-1)*nf+k) += delta2*divU*dNdX(i,k);

      //--- CONTINUITY STABILIZATION ---

      // Momentum residual term
      for (k = 1;k <= nsd;k++) 
	ES(i*nf) -= res(k)*dNdX(i,k);
    }
  }	

  myMats.rhsOnly = false;
  return true;
}
			     

bool NavierStokesG2MP::evalInt (LocalIntegral*& elmInt, double detJW,
				      const Vector& N, const Matrix& dNdX,
				      const Matrix3D& d2NdX2, 
				      const Vec3& X, double h) const
{
  int      i, j;
  int      k, l, m;
  real     diff, divU, divI, divJ, convI, convJ, laplace, U, P;
  real     delta1, delta2; 
  Vector   vel(3), conv(3), res(3), dPdX(3);
  Matrix   dUdX(3,3);
  Matrix3D d2UdX2(3,3,3);

  const int nf = nsd + 1;
  
  Vector& EV = *eVs[0];
  // Velocity in integration point
  vel.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      vel(k) += EV((i-1)*nf + k)*N(i);

  U = 0.0;
  for (k = 1;k <= nsd;k++)
    U += vel(k)*vel(k);
  U = sqrt(U);

  // Pressure in integration point
  P = 0.0;
  for (i = 1;i <= N.size();i++)
    P += EV(i*nf)*N(i);
  
  // Velocity gradient in integration point
  dUdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	dUdX(k,l) += dNdX(i,l)*EV((i-1)*nf+k);
  
  // Pressure gradient in integration point
  dPdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdX(k) += dNdX(i,k)*EV(i*nf);
  
  // Hessian for velocity
  d2UdX2.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	for (m = 1;m <= nsd;m++)
	  d2UdX2(k,m,l) += d2NdX2(i,m,l)*EV((i-1)*nf + k);

  // Stabilization parameters
  if (mu < rho*U*h) {
    // RUNAR
    delta1 = 0.5/sqrt(1.0/(h*h) + (U*U)/(h*h));
    delta2 = h;
  }
  else 
    delta1 = delta2 = h*h;
  delta1 *= detJW;
  delta2 *= detJW;
  // RUNAR
  //delta1 = 0.0001*h*h*detJW;
  //delta2 = 0.0001*h*h*detJW;

  elmInt = &myMats;
  if (eM) {
    Matrix& EM = *eM;

    // Convection terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	// Sum convection for each direction
	convJ = 0.0;
	for (k = 1;k <= nsd;k++)
	  convJ += vel(k)*dNdX(i,k);
	convJ *= rho*N(j)*detJW;
	
	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf+k,(i-1)*nf+k) += convJ;
	
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) += rho*dUdX(l,k)*N(j)*N(i)*detJW;
	}
      }
     
    // Viscous terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= mu*detJW;
	for (k = 1;k <= nsd;k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    if (formulation == STRESS) 
      for (i = 1; i <= N.size(); i++)
 	for (j = 1; j <= N.size(); j++)
 	  for (k = 1; k <= nsd; k++)
 	    for (l = 1; l <= nsd; l++)
 	      EM(nf*(j-1)+l,nf*(i-1)+k) += mu*dNdX(i,k)*dNdX(j,l)*detJW;
    
    // Pressure and divergence terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++)
	for (k = 1;k <= nsd;k++) {
	  divI = N(j)*dNdX(i,k)*detJW;
	  EM(nf*(i-1)+k,nf*j) -= divI;
	  EM(nf*j,nf*(i-1)+k) -= divI;
	}
    
    // Stabilization terms
    
    // Computation of non-linear residual
    res.fill(0.0);
    for (k = 1;k <= nsd;k++) {
      for (l = 1;l <= nsd;l++)
	res(k) += rho*vel(l)*dUdX(k,l) - mu*d2UdX2(k,l,l);

      res(k) += dPdX(k);
      res(k) *= delta1;
    }

    // STABILIZATION TERMS
    for (i = 1;i <= N.size();i++) {
      // Convection for solution basis functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
	convI += vel(k)*dNdX(i,k);
      convI *= rho;
      
      for (j = 1;j <= N.size();j++) {
	// Convection for test functions
	convJ = 0.0;
	for (k = 1;k <= nsd;k++)
	  convJ += vel(k)*dNdX(j,k);
	convJ *= delta1*rho;
	
	//--- MOMENTUM STABILIZATION ---

	// Convection terms
	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf + k,(i-1)*nf + k) += convI*convJ;
	  
	  for (l = 1;l <= nsd;l++) {
	    // RUNAR
	    //EM((j-1)*nf + l,(i-1)*nf + k) += rho*res(l)*dNdX(j,k)*N(i);	 
	    EM((j-1)*nf + l,(i-1)*nf + k) += convJ*rho*dUdX(l,k)*N(i);	 
	  }
	}
	
	// Pressure terms
	for (l = 1;l <= nsd;l++)
	  EM((j-1)*nf + l,i*nf) += dNdX(i,l)*convJ;
	
	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= mu;

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nf + k,(i-1)*nf + k) -= diff*convJ;
	
 	// Continuity-continuity stabilization
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) += delta2*dNdX(i,k)*dNdX(j,l);

	//--- CONTINUITY STABILIZATION ---
	
	// Convection terms
	for (k = 1;k <= nsd;k++)
	  EM(j*nf,(i-1)*nf+k) += delta1*convI*dNdX(j,k);
	
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM(j*nf,(i-1)*nf+k) += delta1*rho*dUdX(l,k)*N(i)*dNdX(j,l);
	
	// Pressure terms
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= delta1;
	EM(j*nf,i*nf) += laplace;
	
	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= delta1*mu;
	
	for (k = 1;k <= nsd;k++)
	  EM(j*nf,(i-1)*nf+k) -= diff*dNdX(j,k);
      }
    }
  }    

  if (eS) {
    // Integrate rhs vector
    Vector& ES = *eS;
    Vector fb(g,nsd);
    fb *= this->getMassDensity(X)*detJW;
    for (i = 1;i <= N.size();i++) {
      // --- Residual of Galerkin formulation ---

      // -- Momentum equations --

      // Body forces
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) += fb(k)*N(i);

      // Convection terms
      for (k = 1;k <= nsd;k++) {
	conv.fill(0.0);
	for (m = 1;m <= nsd;m++)
	  conv(k) += vel(m)*dUdX(k,m);
	
	ES((i-1)*nf + k) -= rho*conv(k)*N(i)*detJW;
      }

      // Pressure terms
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf + k) += P*dNdX(i,k)*detJW;

      // Viscous terms
      for (k = 1;k <= nsd;k++)
	for (m = 1;m <= nsd;m++)
	  ES((i-1)*nf + k) -= mu*dUdX(k,m)*dNdX(i,m)*detJW;

      // -- Continuity equation --
      divU = 0.0;
      for (m = 1;m <= nsd;m++)
	divU += dUdX(m,m);
      
      ES(i*nf) += divU*N(i)*detJW;

      //--- MOMENTUM STABILIZATION---

      // Convection for solution basis functions
      convI = 0.0;
      for (k = 1;k <= nsd;k++)
	convI += vel(k)*dNdX(i,k);
      convI *= rho;

      // Convection terms
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf+k) -= res(k)*convI;

      // Divergence of velocity field and velocity test function
      divU = 0.0;
      for (k = 1;k <= nsd;k++) 
	divU += dUdX(k,k);
      
      // Continuity residual term
      for (k = 1;k <= nsd;k++)
      	ES((i-1)*nf+k) -= delta2*divU*dNdX(i,k);

      //--- CONTINUITY STABILIZATION ---

      // Momentum residual term
      for (k = 1;k <= nsd;k++)
	ES(i*nf) -= res(k)*dNdX(i,k);
    }
  }	

  myMats.rhsOnly = false;
  return true;
}
