// $Id: NavierStokesG2.C,v 1.4 2011-02-08 12:38:23 rho Exp $
//==============================================================================
//!
//! \file NavierStokesG2.C
//!
//! \date May 28 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for G2 stabilized Navier-Stokes problems.
//!
//==============================================================================

#include "NavierStokesG2.h"
#include "TimeDomain.h"
#include "Utilities.h"


NavierStokesG2::NavierStokesG2(short int n, ProblemFormulation form, int itg)
  : StabilizedStokes(n,form,itg)
{
  mu  = 1.0e-3;
  rho = 1.0;

  // Need current solution and solution at previous timestep
  eVs.resize(2);
  for (int i = 0;i < 2;i++)
    eVs[i] = new Vector();
  primsol.resize(2);
}


NavierStokesG2::~NavierStokesG2()
{
}


bool NavierStokesG2::evalInt(LocalIntegral*& elmInt, 
			     const TimeDomain& time, double detJW,
			     const Vector& N, const Matrix& dNdX,
			     const Matrix3D& d2NdX2, const Vec3& X,
			     double h) const
{
  int      i, j;
  int      k, l, m;
  real     diff, divU, divI, divJ, convI, convJ, laplace, U, P;
  real     delta1, delta2; 
  real     dtInv, dUdt;
  Vector   vel(3), pvel(3), conv(3), res(3), dPdX(3);
  Matrix   dUdX(3,3);
  Matrix3D d2UdX2(3,3,3);

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
    if (time.dt == 0.0)
      delta1 = 0.5/sqrt(1.0/(h*h) + (U*U)/(h*h));
    else
      delta1 = 0.5/sqrt(dtInv*dtInv + (U*U)/(h*h));

    delta2 = h;
  }
  else 
    delta1 = delta2 = h*h;
  delta1 *= 0.000001*detJW;
  delta2 *= 0.000001*detJW;

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
			     

bool NavierStokesG2::evalInt (LocalIntegral*& elmInt, double detJW,
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


bool NavierStokesG2::evalSol (Vector& s, const Vector& N,
			      const Matrix& dNdX, const Vec3& X,
			      const std::vector<int>& MNPC) const
{
  if (primsol.empty())
  {
    std::cerr <<" *** NavierStokesG2::evalSol: No primary solution."
              << std::endl;
    return false;
  }

  int ierr = 0;
  if (ierr = utl::gather(MNPC,nsd,primsol[0],*eVs[0]))
  {
    std::cerr <<" *** NavierStokesG2::evalSol: Detected "
              << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

  return true;
}
