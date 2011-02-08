//$Id: ChorinVelPred.C,v 1.3 2011-02-08 12:32:28 rho Exp $
//==============================================================================
//!
//! \file ChorinVelPred.C
//!
//! \date Sep 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for velocity prediction in Chorin's method
//!
//==============================================================================

#include "ChorinVelPred.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Utilities.h"


ChorinVelPred::ChorinVelPred(short int n, ProblemFormulation form, 
			     int itg, bool incPress, bool mixed)
  : Stokes(n,(SIM::Formulation)form,itg), incPressure(incPress), mixedFEM(mixed)
{
  // Number of fields equals number of space dimensions
  nf = nsd;
}


ChorinVelPred::~ChorinVelPred()
{
  // Free memory for pressure element solution vectors
  if (!ePs.empty())
    for (int n = 0;n < ePs.size();n++)
      if (ePs[n]) delete ePs[n];
}


bool ChorinVelPred::initElement(const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eM)  eM->resize(nsd*nen,nsd*nen,true);
  if (eS)  eS->resize(nsd*nen,true);

  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    for (int i = 0;i < primsol.size();i++) 
      if (ierr = utl::gather(MNPC,nsd,primsol[i],*eVs[i]))
        std::cerr <<" *** ChorinVelPred::initElement: Detected "
                  << ierr <<" node numbers out of range."<< std::endl;
    
  if (!ePs.empty() && !psol.empty())
    for (int i = 0;i < psol.size();i++) 
      if (ierr = utl::gather(MNPC,1,psol[i],*ePs[i]))
        std::cerr <<" *** ChorinVelPred::initElement: Detected "
                  << ierr <<" node numbers out of range."<< std::endl;
    
  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinVelPred::initElement(const std::vector<int>& MNPC1,
				const std::vector<int>& MNPC2,
				size_t n1)
{
  // Only velocity degrees of freedom
  const size_t nen1 = MNPC1.size();

  if (eM) eM->resize(nsd*nen1,nsd*nen1,true);
  if (eS) eS->resize(nsd*nen1,true);

  // Extract element velocity vectors
  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    for (int i = 0;i < primsol.size();i++) 
      if (ierr = utl::gather(MNPC1,nsd,primsol[i],*eVs[i]))
        std::cerr <<" *** ChorinVelPredMixed::initElement: Detected "
                  << ierr <<" node numbers out of range."<< std::endl;

  // Extract element pressure vectors
  if (!ePs.empty() && !psol.empty())
    for (int i = 0;i < psol.size();i++)
      if (ierr = utl::gather(MNPC2,1,psol[i],*ePs[i]))
        std::cerr <<" *** ChorinVelPredMixed::initElement: Detected "
                  << ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinVelPred::initElementBou(const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eS) eS->resize(nsd*nen,true);

  myMats->withLHS = false;
  return true;
}


bool ChorinVelPred::initElementBou(const std::vector<int>& MNPC1,
				   const std::vector<int>& MNPC2, 
				   size_t)
{
  // Only velocity degrees of freedom
  const size_t nen1 = MNPC1.size();

  if (eS) eS->resize(nsd*nen1,true);

  myMats->withLHS = false;
  return true;
}



NormBase* ChorinVelPred::getNormIntegrand (AnaSol* asol) const
{
  return new ChorinStokesNorm(*const_cast<ChorinVelPred*>(this),asol);
}


bool ChorinVelPred::evalSol (Vector& s, const Vector& N,
			     const Matrix& dNdX, const Vec3& X,
			     const std::vector<int>& MNPC) const
{
  int ierr = 0;
  if (!eVs.empty() && eVs[0] && !primsol.front().empty())
    if (ierr = utl::gather(MNPC,nsd,primsol.front(),*(eVs[0]))) {
      std::cerr << " *** ChorinVelPred::evalSol: Detected "
		<< ierr << " node numbers out of range." << std::endl;
      return false;
    }
  if (!ePs.empty() && ePs[0] && !psol.front().empty())
    if (ierr = utl::gather(MNPC,1,psol.front(),*(ePs[0]))) {
      std::cerr << " *** ChorinVelPred::evalSol: Detected "
		<< ierr << " node numbers out of range." << std::endl;
      return false;
    } 
  
  // Evaluate stress tensor, sigma  
  SymmTensor sigma(nsd);
  if (!this->stress(N,dNdX,sigma))
    return false;

  s = sigma;
  return true;
}


bool ChorinVelPred::evalSol(Vector& s,
			    const Vector& N1, const Vector& N2,
			    const Matrix& dN1dX, const Matrix& dN2dX, 
			    const Vec3& X,
			    const std::vector<int>& MNPC1,
			    const std::vector<int>& MNPC2) const
{
  // Extract element velocity vector
  int ierr = 0;
  if (!eVs.empty() && eVs[0] && !primsol.front().empty())
    if (ierr = utl::gather(MNPC1,nsd,primsol.front(),*(eVs[0]))) {
      std::cerr << " *** ChorinVelPredMixed::evalSol: Detected "
                << ierr << " node numbers out of range." << std::endl;
      return false;
    }

  // Extract element pressure vector
  if (!ePs.empty() && ePs[0] && !psol.front().empty())
    if (ierr = utl::gather(MNPC2,1,psol.front(),*(ePs[0]))) {
      std::cerr << " *** ChorinVelPredMixed::evalSol: Detected "
                << ierr << " node numbers out of range." << std::endl;
      return false;
    } 

  // Evaluate stress tensor, sigma
  SymmTensor sigma(nsd);
  if (!this->stress(N1,N2,dN1dX,dN2dX,sigma))
    return false;
 
  s = sigma;
  return true;
}


bool ChorinVelPred::evalSol (Vector& s, const TensorFunc& asol, 
                      const Vec3& X) const
{
  s = asol(X);
  return true;
}


bool ChorinVelPred::strain(const Matrix& dNdX, SymmTensor& eps) const
{
  int i, k, l;

  if (dNdX.cols() < nsd) {
    std::cerr <<" *** ChorinVelPred::strain: Invalid dimension on dNdX "
              << dNdX.rows() <<" "<< dNdX.cols() << std::endl;
    return false;
  }

  if (!eVs.empty() && !eVs[0]->empty() && eps.dim() > 0) {
    Vector& EV = *(eVs[0]);
    for (i = 1;i <= dNdX.rows();i++)
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= k;l++)
	  eps(k,l) = dNdX(i,l)*EV((i-1)*nsd+k); + dNdX(i,k)*EV((i-1)*nsd+l);
    
    eps *= mu;
  }

  return true;
}


bool ChorinVelPred::stress(const Vector& N, const Matrix& dNdX, SymmTensor& sigma) const
{
  int i, k;

  // Compute strain
  if (!this->strain(dNdX,sigma))
    return false;
  
  // Compute pressure
  if (!ePs.empty() && !ePs[0]->empty() && sigma.dim() > 0) {
    double P = 0.0;
    Vector& EP = *(ePs[0]);
    for (i = 1;i <= N.size();i++)
      P += EP(i)*N(i);
    
    // Add pressure and strain contributions
    for (k = 1;k <= nsd;k++)
      sigma(k,k) -= P;
  }

  return true;
}


bool ChorinVelPred::stress(const Vector& N1, const Vector& N2,
			   const Matrix& dN1dX, const Matrix& dN2dX,
			   SymmTensor& sigma) const
{
  int i, k;

  // Compute strain
  if (!this->strain(dN1dX,sigma))
    return false;

  // Add pressure
  if (!ePs.empty() && ePs[0]) {
    double P = 0.0;
    Vector& EP = *ePs[0];
    for (i = 1;i <= N2.size();i++)
      P += EP(i)*N2(i);

    // Add pressure to strain
    for (k = 1;k <= nsd;k++)
      sigma(k,k) -= P;
  }
                
  return true;
}


bool ChorinStokesNorm::evalInt (LocalIntegral*& elmInt, double detJW,
				const Vector& N, const Matrix& dNdX,
				const Vec3& X) const
{
  int    i, k, l;
  double value, eps, epsh;

  const int nsd = dNdX.cols();
  const int nf  = nsd+1;
  const int nbf = N.size();

  if (!anasol->hasScalarSol() && !anasol->hasVectorSol()) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No analytical solution."<< std::endl;
    return false;
  }

  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (!eNorm) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No norm array."<< std::endl;
    return false;
  }

  // Find current element velocity vector
  Vector* eV = const_cast<Vector*>(problem.getElementVelocity());
  if (!eV) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No velocity" << std::endl;
    return false;
  }
  Vector& EV = *eV;

  // Find current element pressure vector
  Vector* eP = const_cast<Vector*>(problem.getElementPressure());
  if (!eP) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No pressure" << std::endl;
    return false;
  }
  Vector& EP = *eP;

  ElmNorm& pnorm = *eNorm;
  int ip = 0;
  // Pressure norms
  if (anasol->getScalarSol()) {
    // Analytical pressure
    real P = (*anasol->getScalarSol())(X);

    // Computed pressure
    real Ph = 0.0;
    for (i = 1;i <= nbf;i++)
      Ph += EP(i)*N(i);

    // Integral of computed pressure 
    pnorm[ip++] += Ph*detJW;
    
    // Integral of basis functions
    pnorm[ip++] += detJW;

    // L2-norm of pressure error
    P -= Ph;
    pnorm[ip++] += P*P*detJW;
  }

  // Velocity norms
  if (anasol->getVectorSol()) {
    // Analytical velocity
    Vec3 U = (*anasol->getVectorSol())(X);

    // Computed velocity
    Vec3 Uh; Uh = 0.0;
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	Uh[k-1] += EV((i-1)*nsd + k)*N(i);
    
    // L2-norm of velocity error
    U -= Uh;
    pnorm[ip++] += U*U*detJW;
  }

  // Pressure gradient norm
  if (anasol->getScalarSecSol()) {
    // Analytical pressure gradient
    Vec3 dP = (*anasol->getScalarSecSol())(X);
    
    // Computed pressure gradient
    Vec3 dPh; dPh = 0.0;
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	dPh[k-1] += EP(i)*dNdX(i,k);

    // H1-seminorm of pressure
    dP -= dPh;
    pnorm[ip++] += dP*dP*detJW;
  }

  // Energy norms
  if (anasol->getVectorSecSol()) {
    // Viscosity
    const double mu = problem.getViscosity(X);

    // Analytical velocity gradient
    Tensor gradU = (*anasol->getVectorSecSol())(X);
    
    // Computed velocity gradient
    Matrix gradUh(nsd,nsd);
    gradUh.fill(0.0);
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++)
	  gradUh(k,l) += EV((i-1)*nsd + k)*dNdX(i,l);

    if (problem.getFormulation() == Stokes::LAPLACE) {
      for (k = 1;k <= nsd;k++) 
	for (l = 1;l <= nsd;l++) {
	  // Energy norm of analytical solution
	  pnorm[ip] += mu*gradU(k,l)*gradU(k,l)*detJW;
	  // Energy norm of numerical solution
	  pnorm[ip+1] += mu*gradUh(k,l)*gradUh(k,l)*detJW;
	  // Energy norm of error
	  value = gradU(k,l)-gradUh(k,l);
	  pnorm[ip+2] += mu*value*value*detJW;
	}
    }
    else {
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++) {
	  // Strain of analytical solution
	  eps  = 0.5*(gradU(k,l) + gradU(l,k));
	  // Strain of computed solution
	  epsh = 0.5*(gradUh(k,l) + gradUh(l,k));
	  // Energy norm of analytical solution
	  pnorm[ip] += mu*eps*eps*detJW;
	  // Energy norm of computed solution
	  pnorm[ip+1] += mu*epsh*epsh*detJW;
	  // Energy norm of error
	  eps -= epsh;
	  pnorm[ip+2] += mu*epsh*epsh*detJW;
	} 
    }
  }

  return true;
}


bool ChorinStokesNorm::evalInt(LocalIntegral*& elmInt,
			       const TimeDomain& time, double detJW,
			       const Vector& N1, const Vector& N2,
			       const Matrix& dN1dX, const Matrix& dN2dX,
			       const Vec3& X) const
{
  int    i, k, l;
  double value, eps, epsh;

  const int nsd  = dN1dX.cols();
  const int nf   = nsd+1;
  const int nbfU = N1.size();
  const int nbfP = N2.size();

  if (!anasol->hasScalarSol() && !anasol->hasVectorSol()) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No analytical solution."<< std::endl;
    return false;
  }

  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (!eNorm) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No norm array."<< std::endl;
    return false;
  }

  // Find current element velocity vector
  Vector* eV = const_cast<Vector*>(problem.getElementVelocity());
  if (!eV) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No svelocity" << std::endl;
    return false;
  }
  Vector& EV = *eV;

  // Find current element pressure vector
  Vector* eP = const_cast<Vector*>(problem.getElementPressure());
  if (!eP) {
    std::cerr <<" *** ChorinStokesNorm::evalInt: No pressure" << std::endl;
    return false;
  }
  Vector& EP = *eP;

  ElmNorm& pnorm = *eNorm;
  int ip = 0;
  // Pressure norms
  if (anasol->getScalarSol()) {
    // Analytical pressure
    real P = (*anasol->getScalarSol())(X);

    // Computed pressure
    real Ph = 0.0;
    for (i = 1;i <= nbfP;i++)
      Ph += EP(i)*N2(i);

    // Integral of computed pressure 
    pnorm[ip++] += Ph*detJW;
    
    // Integral of basis functions
    pnorm[ip++] += detJW;

    // L2-norm of pressure error
    P -= Ph;
    pnorm[ip++] += P*P*detJW;
  }

  // Velocity norms
  if (anasol->getVectorSol()) {
    // Analytical velocity
    Vec3 U = (*anasol->getVectorSol())(X);

    // Computed velocity
    Vec3 Uh; Uh = 0.0;
    for (i = 1;i <= nbfU;i++)
      for (k = 1;k <= nsd;k++)
	Uh[k-1] += EV((i-1)*nsd + k)*N1(i);
    
    // L2-norm of velocity error
    U -= Uh;
    pnorm[ip++] += U*U*detJW;
  }

  // Pressure gradient norm
  if (anasol->getScalarSecSol()) {
    // Analytical pressure gradient
    Vec3 dP = (*anasol->getScalarSecSol())(X);
    
    // Computed pressure gradient
    Vec3 dPh; dPh = 0.0;
    for (i = 1;i <= nbfP;i++)
      for (k = 1;k <= nsd;k++)
	dPh[k-1] += EP(i)*dN2dX(i,k);

    // H1-seminorm of pressure
    dP -= dPh;
    pnorm[ip++] += dP*dP*detJW;
  }

  // Energy norms
  if (anasol->getVectorSecSol()) {
    // Viscosity
    const double mu = problem.getViscosity(X);

    // Analytical velocity gradient
    Tensor gradU = (*anasol->getVectorSecSol())(X);
    
    // Computed velocity gradient
    Matrix gradUh(nsd,nsd);
    gradUh.fill(0.0);
    for (i = 1;i <= nbfU;i++)
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++)
	  gradUh(k,l) += EV((i-1)*nsd + k)*dN1dX(i,l);

    if (problem.getFormulation() == Stokes::LAPLACE) {
      for (k = 1;k <= nsd;k++) 
	for (l = 1;l <= nsd;l++) {
	  // Energy norm of analytical solution
	  pnorm[ip] += mu*gradU(k,l)*gradU(k,l)*detJW;
	  // Energy norm of numerical solution
	  pnorm[ip+1] += mu*gradUh(k,l)*gradUh(k,l)*detJW;
	  // Energy norm of error
	  value = gradU(k,l)-gradUh(k,l);
	  pnorm[ip+2] += mu*value*value*detJW;
	}
    }
    else {
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++) {
	  // Strain of analytical solution
	  eps  = 0.5*(gradU(k,l) + gradU(l,k));
	  // Strain of computed solution
	  epsh = 0.5*(gradUh(k,l) + gradUh(l,k));
	  // Energy norm of analytical solution
	  pnorm[ip] += mu*eps*eps*detJW;
	  // Energy norm of computed solution
	  pnorm[ip+1] += mu*epsh*epsh*detJW;
	  // Energy norm of error
	  eps -= epsh;
	  pnorm[ip+2] += mu*epsh*epsh*detJW;
	} 
    }
  }

  return true;
}


