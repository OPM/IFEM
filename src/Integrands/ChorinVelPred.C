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

  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    if (ierr = utl::gather(MNPC,nsd,primsol[0],*eVs[0]))
      std::cerr <<" *** ChorinVelPred::initElementBou: Detected "
                << ierr <<" node numbers out of range."<< std::endl;
    
  if (!ePs.empty() && !psol.empty())
    if (ierr = utl::gather(MNPC,1,psol[0],*ePs[0]))
      std::cerr <<" *** ChorinVelPred::initElementBou: Detected "
                << ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


bool ChorinVelPred::initElementBou(const std::vector<int>& MNPC1,
				   const std::vector<int>& MNPC2, 
				   size_t)
{
  // Only velocity degrees of freedom
  const size_t nen1 = MNPC1.size();

  if (eS) eS->resize(nsd*nen1,true);

  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    if (ierr = utl::gather(MNPC1,nsd,primsol[0],*eVs[0]))
      std::cerr <<" *** ChorinVelPred::initElementBou: Detected "
                << ierr <<" node numbers out of range."<< std::endl;
    
  if (!ePs.empty() && !psol.empty())
    if (ierr = utl::gather(MNPC2,1,psol[0],*ePs[0]))
      std::cerr <<" *** ChorinVelPred::initElementBou: Detected "
                << ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


NormBase* ChorinVelPred::getNormIntegrand (AnaSol* asol) const
{
  return new ChorinStokesNorm(*const_cast<ChorinVelPred*>(this),asol);
}


NormBase* ChorinVelPred::getForceIntegrand (AnaSol* asol) const
{
  return new ChorinStokesForce(*const_cast<ChorinVelPred*>(this),asol);
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


bool ChorinVelPred::strain(const Matrix& dNdX, Tensor& eps) const
{
  int i, k, l;

  if (dNdX.cols() < nsd) {
    std::cerr <<" *** ChorinVelPred::strain: Invalid dimension on dNdX "
              << dNdX.rows() <<" "<< dNdX.cols() << std::endl;
    return false;
  }

  eps.zero();
  if (!eVs.empty() && !eVs[0]->empty() && eps.dim() > 0) {
    Vector& EV = *(eVs[0]);
    for (i = 1;i <= dNdX.rows();i++)
      for (k = 1;k <= nsd;k++)
        for (l = 1;l <= nsd;l++) 
          eps(k,l) += dNdX(i,l)*EV((i-1)*nsd+k);
  }
  else 
    return false;
    
  if (formulation == Stokes::STRESS) {
    eps += eps.transpose();
    eps *= 0.5;
  }

  return true;
}


bool ChorinVelPred::stress(const Vector& N, const Matrix& dNdX, Tensor& sigma) const
{
  int i, k;

  // Compute strain
  if (!this->strain(dNdX,sigma))
    return false;
  sigma *= mu;
  
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
			   Tensor& sigma) const
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
    
    // Numerical velocity gradient
    Tensor gradUh(3);
    gradUh.zero();
    problem.strain(dNdX,gradUh);

    if (problem.getFormulation() == Stokes::STRESS) { 
      pnorm[ip++] += 2.0*mu*gradU.innerProd(gradU)*detJW;
      pnorm[ip++] += 2.0*mu*gradUh.innerProd(gradUh)*detJW;
      gradUh *= -1.0;
      gradU += gradUh;
      pnorm[ip++] += 2.0*mu*gradU.innerProd(gradU)*detJW; 
    }
    else {
      pnorm[ip++] += mu*gradU.innerProd(gradU)*detJW;
      pnorm[ip++] += mu*gradUh.innerProd(gradUh)*detJW;
      gradUh *= -1.0;
      gradU += gradUh;
      pnorm[ip++] += mu*gradU.innerProd(gradU)*detJW; 
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
    Tensor gradUh(nsd);
    problem.strain(dN1dX,gradUh);
    
    if (problem.getFormulation() == Stokes::STRESS) { 
      pnorm[ip++] += 2.0*mu*gradU.innerProd(gradU)*detJW;
      pnorm[ip++] += 2.0*mu*gradUh.innerProd(gradUh)*detJW;
      gradUh *= -1.0;
      gradU += gradUh;
      pnorm[ip++] += 2.0*mu*gradU.innerProd(gradU)*detJW; 
    }
    else {
      pnorm[ip++] += mu*gradU.innerProd(gradU)*detJW;
      pnorm[ip++] += mu*gradUh.innerProd(gradUh)*detJW;
      gradUh *= -1.0;
      gradU += gradUh;
      pnorm[ip++] += mu*gradU.innerProd(gradU)*detJW; 
    }
  }

  return true;
}


bool ChorinStokesForce::initElementBou(const std::vector<int>& MNPC)
{
  return problem.initElementBou(MNPC);
}


bool ChorinStokesForce::initElementBou(const std::vector<int>& MNPC1,
                                       const std::vector<int>& MNPC2, 
                                       size_t n1)
{
  return problem.initElementBou(MNPC1,MNPC2,n1);
}


bool ChorinStokesForce::evalBou(LocalIntegral*& elmInt,
                          const TimeDomain& time, double detJW,
                          const Vector& N, const Matrix& dNdX,
                          const Vec3& X, const Vec3& normal) const
{
  const int nsd = dNdX.cols();

  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (!eNorm) {
    std::cerr <<" *** StokesForce::evalBou: No force array."<< std::endl;
    return false;
  }  

  ChorinVelPred* prob = static_cast<ChorinVelPred*>(&problem);
  
  double mu = prob->getViscosity(X);

  // Numerical approximation of stress
  Tensor sigmah(3);
  prob->stress(N,dNdX,sigmah);

  // Traction
  Vec3 th = sigmah*normal;

  // Numerical force term
  ElmNorm& pnorm = *eNorm;
  int i, ip = 0;
  for (i = 0;i < nsd;i++)
    pnorm[ip++] += th[i]*detJW;

    // Analytical force term and error norm
  if ( anasol && anasol->getScalarSol() && anasol->getVectorSecSol()) {
    real P = (*anasol->getScalarSol())(X);
    Tensor sigma = (*anasol->getVectorSecSol())(X);

    // Symmetrice for stress formulation
    if (prob->getFormulation() == Stokes::STRESS) 
      sigma += sigma.transpose();

    // Analytical stress
    sigma *= mu;
    for (i = 1;i <= nsd;i++)
      sigma(i,i) -= P;

    // Analytical traction
    Vec3 t = sigma*normal;

    for (i = 0;i < nsd;i++)
      pnorm[ip++] += t[i]*detJW;

    // Error in traction
    t -= th;
    pnorm[ip++] += t*t*detJW;
  }

  return true;
}


