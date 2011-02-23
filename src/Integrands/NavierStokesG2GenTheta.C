//$Id: NavierStokesG2GenTheta.C,v 1.3 2011-02-08 12:40:31 rho Exp $
//==============================================================================
//!
//! \file NavierStokesG2GenTheta.C
//!
//! \date August 3 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for G2 stabilized Navier-Stokes problems.
//! \details Time integration by the generalized theta method.
//!
//==============================================================================

#include "NavierStokesG2GenTheta.h"
#include "TimeDomain.h"
#include "ElmNorm.h"


NavierStokesG2GenTheta::NavierStokesG2GenTheta(short int n, SIM::Formulation form, int itg)
  : NavierStokesG2(n,form,itg) {}


NavierStokesG2GenTheta::~NavierStokesG2GenTheta() {}


bool NavierStokesG2GenTheta::evalInt(LocalIntegral*& elmInt, 
				     const TimeDomain& time, double detJW,
				     const Vector& N, const Matrix& dNdX,
				     const Matrix3D& d2NdX2, const Vec3& X,
				     double h) const
{
  int      i, j;
  int      k, l, m;
  real     temp1, temp2;
  real     diff, divU, divI, divJ, convI, convJ, laplace, U, P;
  real     delta1, delta2; 
  real     dtInv, dUdt;
  Vector   vel(3), pvel(3), conv(3), convp(3), diffp(3), res(3), resp(3), dPdX(3);
  Matrix   dUdX(3,3), dUdXn(3,3);
  Matrix3D d2UdX2(3,3,3), d2UdX2n(3,3,3);

  // Time integration parameters
  const real thetau = thetat;
  const real thetap = thetam;

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
  
  // Velocity gradient at previous timestep
  dUdXn.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	dUdXn(k,l) += dNdX(i,l)*EVP((i-1)*nf+k);

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
  
  // Hessian for velocity at previous timestep
  d2UdX2n.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	for (m = 1;m <= nsd;m++)
	  d2UdX2n(k,m,l) += d2NdX2(i,m,l)*EVP((i-1)*nf + k);

  // Stabilization parameters
  if (mu < rho*U*h) {
    if (time.dt == 0.0)
      delta1 = 0.5/sqrt(1.0/(h*h) + (U*U)/(h*h));
    else
      delta1 = 0.5/sqrt(dtInv*dtInv + (U*U)/(h*h));

    delta2 = h;
  }
  else 
    delta1 = delta2 = h*h;
  delta1 *= 0.01*detJW;
  delta2 *= 0.01*detJW;
  // RUNAR
  //delta1 = 0.0001*h*h*detJW;
  //delta2 = 0.0001*h*h*detJW;

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
	  convJ += vel(k)*dNdX(i,k);
	convJ *= thetau*rho*N(j)*detJW;

	temp1 = thetau*rho*N(i)*N(j)*detJW;
	
	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf+k,(i-1)*nf+k) += convJ;
	
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) += temp1*dUdX(l,k);
	}
      }
     
    // Viscous terms
    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	
	laplace *= thetau*mu*detJW;
	for (k = 1;k <= nsd;k++) 
	  EM(nf*(j-1)+k,nf*(i-1)+k) += laplace;       
      }
    
    temp1 = thetau*mu*detJW;
    if (formulation == SIM::STRESS) 
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
	  EM(nf*(i-1)+k,nf*j) -= thetap*divI;
	  EM(nf*j,nf*(i-1)+k) -= divI;
	}
    
    // Stabilization terms
    
    // Computation of non-linear residual
    res.fill(0.0);
    resp.fill(0.0);
    for (k = 1;k <= nsd;k++) {
      for (l = 1;l <= nsd;l++) {
	res(k)  += rho*vel(l)*dUdX(k,l)   - mu*d2UdX2(k,l,l);
	resp(k) += rho*pvel(l)*dUdXn(k,l) - mu*d2UdX2n(k,l,l);
      }

      //res(k) += dPdX(k);
      res(k)  *= delta1;
      resp(k) *= delta1;
    }

    // Convection at previous timestep
    convp.fill(0.0);
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	convp(k) += rho*pvel(l)*dUdXn(k,l);

    // Viscous terms at previous timestep
    diffp.fill(0.0);
    for (k = 1;k <= nsd;k++)
      for (l = 1;l <= nsd;l++)
	diffp(k) += mu*d2UdX2(k,l,l);

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
	  EM((j-1)*nf + k,(i-1)*nf + k) += thetau*convI*convJ;
	  
	  for (l = 1;l <= nsd;l++) {
	    EM((j-1)*nf + l,(i-1)*nf + k) += thetau*rho*res(l)*dNdX(j,k)*N(i);	 
	    EM((j-1)*nf + l,(i-1)*nf + k) += thetau*convJ*rho*dUdX(l,k)*N(i);	 
	    EM((j-1)*nf + l,(i-1)*nf + k) += delta1*alpha*thetap*convp(l)*rho*N(i)*dNdX(j,k);
	  }
	}
	
	// Pressure terms
	for (l = 1;l <= nsd;l++) {
	  EM((j-1)*nf + l,i*nf) += thetap*dNdX(i,l)*convJ;
	
	  for (k = 1;k <= nsd;k++)
	    EM((j-1)*nf + l,(i-1)*nf + k) += delta1*thetap*dPdX(l)*rho*N(i)*dNdX(j,k); 
	}

	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= thetau*mu;

	for (k = 1;k <= nsd;k++) {
	  EM((j-1)*nf + k,(i-1)*nf + k) -= diff*convJ;

	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf + l,(i-1)*nf + k) -= delta1*alpha*thetap*diffp(l)*rho*N(i)*dNdX(j,k);
	}	

 	// Continuity-continuity stabilization
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM((j-1)*nf+l,(i-1)*nf+k) += delta2*dNdX(i,k)*dNdX(j,l);

	//--- CONTINUITY STABILIZATION ---
	
	// Convection terms
	for (k = 1;k <= nsd;k++)
	  EM(j*nf,(i-1)*nf+k) += thetau*delta1*convI*dNdX(j,k);
	
	for (k = 1;k <= nsd;k++)
	  for (l = 1;l <= nsd;l++)
	    EM(j*nf,(i-1)*nf+k) += thetau*delta1*rho*dUdX(l,k)*N(i)*dNdX(j,l);
	
	// Pressure terms
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= thetap*delta1;
	EM(j*nf,i*nf) += laplace;
	
	// Viscous terms
	diff = 0.0;
	for (k = 1;k <= nsd;k++)
	  diff += d2NdX2(i,k,k);
	diff *= thetau*delta1*mu;
	
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
	ES((i-1)*nf+k) += thetap*fb(k)*N(i);

      // Convection terms
      for (k = 1;k <= nsd;k++) {
	conv.fill(0.0);
	for (m = 1;m <= nsd;m++)
	  conv(k) += vel(m)*dUdX(k,m);
	
	ES((i-1)*nf + k) -= thetau*rho*conv(k)*N(i)*detJW;
      }

      // Pressure terms
      for (k = 1;k <= nsd;k++)
	ES((i-1)*nf + k) += thetap*P*dNdX(i,k)*detJW;

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
      for (k = 1;k <= nsd;k++) {
	ES((i-1)*nf+k) -= thetau*res(k)*convI; 
	ES((i-1)*nf+k) -= alpha*thetap*convp(k)*N(i)*detJW;
	ES((i-1)*nf+k) -= delta1*alpha*thetap*convp(k)*convI;
      }

      // Viscous terms
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++) {
	  ES((i-1)*nf+k) += alpha*thetap*mu*dUdXn(k,l)*dNdX(j,l)*detJW;
	  ES((i-1)*nf+k) += delta1*alpha*thetap*mu*diffp(k)*convI;
	}

      // Divergence of velocity field and velocity test function
      divU = 0.0;
      for (k = 1;k <= nsd;k++) 
	divU += dUdX(k,k);
      
      // Continuity residual term
      for (k = 1;k <= nsd;k++)
      	ES((i-1)*nf+k) -= delta2*divU*dNdX(i,k);

      //--- CONTINUITY STABILIZATION ---

      // Momentum residual term
      for (k = 1;k <= nsd;k++) {
	ES(i*nf) -= thetau*res(k)*dNdX(i,k);
	ES(i*nf) -= alpha*thetap*resp(k)*dNdX(i,k);
      }
    }
  }	

  return getIntegralResult(elmInt);
}
			     

bool NavierStokesG2GenTheta::evalInt (LocalIntegral*& elmInt, double detJW,
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
    
    if (formulation == SIM::STRESS) 
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
	    EM((j-1)*nf + l,(i-1)*nf + k) += rho*res(l)*dNdX(j,k)*N(i);	 
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

  return getIntegralResult(elmInt);
}
