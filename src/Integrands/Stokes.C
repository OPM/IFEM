// $Id: Stokes.C,v 1.3 2011-02-08 12:23:54 rho Exp $
//==============================================================================
//!
//! \file Stokes.C
//!
//! \date Sep 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for Stokes problems.
//!
//==============================================================================

#include "Stokes.h"
#include "Utilities.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "VTF.h"


Stokes::Stokes (unsigned short int n, SIM::Formulation form, int itg_type)
{
  nsd = n;
  nf = n+1;
  formulation = form;
  integrandType = itg_type;

  mu = rho = 1.0;

  tracFld = 0;

  myMats = new ElmMats();
  myMats->resize(1,1);
  myMats->rhsOnly = false;

  eM = &myMats->A[0];
  eS = &myMats->b[0];
}


Stokes::~Stokes ()
{
  if (myMats) delete myMats;

  for (int i = 0; i < 2; i++)
    if (eVs[i]) delete eVs[i];
}


bool Stokes::initElement (const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  eM->resize(nf*nen,nf*nen,true);
  eS->resize(nf*nen,true);

  int ierr = 0;
  for (size_t i = 0; i < 2 && i < primsol.size() && eVs[i]; i++)
    if (ierr = utl::gather(MNPC,nf,primsol[i],*eVs[i]))
      std::cerr <<" *** Stokes::initElement: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool Stokes::initElementBou (const std::vector<int>& MNPC)
{
  eS->resize(nf*MNPC.size(),true);

  myMats->withLHS = false;
  return true;
}


bool Stokes::evalBou (LocalIntegral*& elmInt, double detJW,
		      const Vector& N, const Matrix& dNdX,
		      const Vec3& X, const Vec3& normal) const
{
  if (!tracFld)
  {
    std::cerr <<" *** Stokes::evalBou: No tractions."<< std::endl;
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
      ES(nf*(a-1)+d) += T[d-1]*N(a)*detJW;

  return getIntegralResult(elmInt);
}


bool Stokes::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
{
  if (tracVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write boundary tractions as discrete point vectors to the VTF-file
  if (!vtf->writeVectors(tracVal,++nBlock))
    return false;

  return vtf->writeVblk(nBlock,"Tractions",1,iStep);
}


bool Stokes::evalSol (Vector& s, const Vector& N,
		      const Matrix& dNdX, const Vec3& X,
		      const std::vector<int>& MNPC) const
{
  return false;
}


bool Stokes::evalSol (Vector& s, const TensorFunc& asol, const Vec3& X) const
{
  s = asol(X);
  return true;
}


const char* Stokes::getFieldLabel(size_t i, const char* prefix) const
{
  static const char* s[3] = { "u_x","u_y","u_z" };

  static std::string label;
  if (prefix)
    label = prefix + std::string(" ");
  else
    label.clear();

  label += s[i];
  return label.c_str();
}


NormBase* Stokes::getNormIntegrand (AnaSol* asol) const
{
  return new StokesNorm(*const_cast<Stokes*>(this),asol);
}


bool Stokes::getIntegralResult (LocalIntegral*& elmInt) const
{
  elmInt = myMats;
  return elmInt ? true : false;
}


bool StokesNorm::initElement (const std::vector<int>& MNPC)
{
  return problem.initElement(MNPC);
}


bool StokesNorm::evalInt (LocalIntegral*& elmInt, double detJW,
			  const Vector& N, const Matrix& dNdX,
			  const Vec3& X) const
{
  int    i, k, l;
  double value, eps, epsh;

  const int nsd = dNdX.cols();
  const int nf  = nsd+1;
  const int nbf = N.size();

  if (!anasol->hasScalarSol() && !anasol->hasVectorSol()) {
    std::cerr <<" *** StokesNorm::evalInt: No analytical solution."<< std::endl;
    return false;
  }

  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (!eNorm) {
    std::cerr <<" *** StokesNorm::evalInt: No norm array."<< std::endl;
    return false;
  }

  // Find current element solution vector
  Vector* eV = const_cast<Vector*>(problem.getElementSolution());
  if (!eV) {
    std::cerr <<" *** StokesNorm::evalInt: No solution"<< std::endl;
    return false;
  }
  Vector& EV = *eV;

  ElmNorm& pnorm = *eNorm;
  int ip = 0;
  // Pressure norms
  if (anasol->getScalarSol()) {
    // Analytical pressure
    real P = (*anasol->getScalarSol())(X);

    // Computed pressure
    real Ph = 0.0;
    for (i = 1;i <= nbf;i++)
      Ph += EV(i*nf)*N(i);

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
	Uh[k-1] += EV((i-1)*nf + k)*N(i);
    
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
	dPh[k-1] += EV(i*nf)*dNdX(i,k);

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
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	for (l = 1;l <= nsd;l++)
	  gradUh(k,l) += EV((i-1)*nf + k)*dNdX(i,l);
    
    if (problem.getFormulation() == Stokes::LAPLACE) {
      for (k = 1;k <= nsd;k++) 
	for (l = 1;l <= nsd;l++) {
	  // Energy norm of analytical solution
	  pnorm[ip++] += mu*gradU(k,l)*gradU(k,l)*detJW;
	  // Energy norm of numerical solution
	  pnorm[ip++] += mu*gradUh(k,l)*gradUh(k,l)*detJW;
	  // Energy norm of error
	  value = gradU(k,l)-gradUh(k,l);
	  pnorm[ip++] += mu*value*value*detJW;
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
	  pnorm[ip++] += mu*eps*eps*detJW;
	  // Energy norm of computed solution
	  pnorm[ip++] += mu*epsh*epsh*detJW;
	  // Energy norm of error
	  eps -= epsh;
	  pnorm[ip++] += mu*epsh*epsh*detJW;
	} 
    }
  }

  return true;
}


bool Stokes::getBodyForce(const Vec3& X, Vector& f) const
{
  const Vec4* Y = dynamic_cast<const Vec4*>(&X);
  if (!Y) return false;

  const double PI = 3.141592653589793238462;
  const double x  = X[0];
  const double y  = X[1];
  const double t  = Y->t;

  f.fill(0.0);
  

  f(1) = rho*pow(sin(PI*x),2.0)*sin(2.0*PI*y)*cos(t) + 
    rho*PI*pow(sin(PI*x)*sin(2.0*PI*y)*sin(t),2.0)*sin(2.0*PI*x) -
    2.0*rho*PI*pow(sin(PI*x)*sin(PI*y)*sin(t),2.0)*sin(2.0*PI*x)*cos(2.0*PI*y) -
    PI*sin(PI*x)*sin(PI*y)*sin(t) - 
    2.0*pow(PI,2.0)*mu*sin(2.0*PI*y)*sin(t)*(cos(2*PI*x)-2*pow(sin(PI*x),2.0));

  f(2) = -rho*sin(2.0*PI*x)*pow(sin(PI*y),2.0)*cos(t) -
    2.0*rho*PI*pow(sin(PI*x)*sin(PI*y)*sin(t),2.0)*cos(2.0*PI*x)*sin(2.0*PI*y) +
    rho*PI*pow(sin(2.0*PI*x)*sin(PI*y)*sin(t),2.0)*sin(2.0*PI*y) +
    PI*cos(PI*x)*cos(PI*y)*sin(t) -
    4*mu*pow(PI*sin(PI*y),2.0)*sin(2.0*PI*x)*sin(t) +
    2.0*mu*pow(PI,2.0)*sin(2.0*PI*x)*cos(2.0*PI*y)*sin(t);

  return true;
}
