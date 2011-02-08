//==============================================================================
//!
//! \file ChorinPressCorr.C
//!
//! \date Sep 30 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for pressure correction in Chorin's method
//!
//==============================================================================

#include "ChorinPressCorr.h"
#include "Utilities.h"
#include "Vec3Oper.h"

ChorinPressCorr::ChorinPressCorr(unsigned short int n, double B0, 
				 bool incPress, bool mixed) 
  : nsd(n), Beta0(B0), incPressure(incPress), mixedFEM(mixed)
{
  // Set default values for density and viscosity
  rho = mu = 1.0;

  // Initialize pointers to zero
  tracFld = 0;

  // Initialize element matrix and vector
  myMats = new ElmMats();
  myMats->A.resize(1);
  myMats->b.resize(1);
  myMats->rhsOnly = false;

  eM = &myMats->A[0];
  eS = &myMats->b[0];

  // Current velocity and pressure solutions are needed
  int nUsols = 1;
  usol.resize(nUsols);
  eVs.resize(nUsols);
  for (int i = 0;i < nUsols;i++)   
    eVs[i] = new Vector();
  int nPsols = 1;
  primsol.resize(nPsols);
  ePs.resize(nPsols);
  for (int i = 0;i < nPsols;i++)   
    ePs[i] = new Vector();
}


ChorinPressCorr::~ChorinPressCorr()
{
  int i;

  // Delete element matrix/vector
  if (myMats) delete myMats;

  // Delete element solution vectors
  if (!eVs.empty())
    for (i = 0;i < eVs.size();i++)
      if (eVs[i]) delete eVs[i];

  if (!ePs.empty())
    for (i = 0;i < ePs.size();i++)
      if (ePs[i]) delete ePs[i];
}


bool ChorinPressCorr::initElement(const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eM) eM->resize(nen,nen,true);
  if (eS) eS->resize(nen,true);

  int ierr = 0;
  if (!ePs.empty() && !primsol.empty())
    for (int i = 0;i < ePs.size();i++)
      if (ierr = utl::gather(MNPC,1,primsol[i],*ePs[i]))
	std::cerr << "*** ChorinPressCorr::initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;
  if (!eVs.empty() && !usol.empty())
    for (int i = 0;i < eVs.size();i++)
      if (ierr = utl::gather(MNPC,nsd,usol[i],*eVs[i]))
	std::cerr << "*** ChorinPressCorr::initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinPressCorr::initElement(const std::vector<int>& MNPC1,
				  const std::vector<int>& MNPC2,
				  size_t n1)
{
  // Only pressure degrees of freedom
  const size_t nen2 = MNPC2.size();

  if (eM) eM->resize(nen2,nen2,true);
  if (eS) eS->resize(nen2,true);

  // RUNAR
  std::vector<int> MNPC = MNPC2;
  for (int i = 0;i < MNPC.size();i++)
    MNPC[i] -= n1;

  // Extract element pressure vectors
  int ierr = 0;
  if (!ePs.empty() && !primsol.empty())
    for (int i = 0;i < ePs.size();i++)
      if (ierr = utl::gather(MNPC2,1,primsol[i],*ePs[i]))
	std::cerr << "*** ChorinPressCorrMixed:initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;

  // Extract element velocity vectors
  if (!eVs.empty() && !usol.empty()) 
    for (int i = 0;i < eVs.size();i++) 
      if (ierr = utl::gather(MNPC1,nsd,usol[i],*eVs[i])) 
	std::cerr << "*** ChorinPressCorrMixed::initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinPressCorr::initElementBou(const std::vector<int>& MNPC)
{				     
  const size_t nen = MNPC.size();
  
  if (eS) eS->resize(nsd*nen,true);
  
  myMats->withLHS = false;
  return true;
}


bool ChorinPressCorr::initElementBou(const std::vector<int>& MNPC1,
				     const std::vector<int>& MNPC2, 
				     size_t n1)
{
  // Only pressure degrees of freedom
  const size_t nen2 = MNPC2.size();

  if (eS) eS->resize(nen2,true);

  // Extract the element level velocity vector
  // Extract element pressure vectors
  int ierr = 0;
  if (!ePs.empty() && !primsol.empty())
    if (ierr = utl::gather(MNPC2,1,primsol[0],*ePs[0]))
      std::cerr << "*** ChorinPressCorrMixed:initElement: Detected "
		<< ierr << " node numbers out of range." << std::endl;

  // Extract element velocity vectors
  if (eVs.empty() && !usol.empty())
    if (ierr = utl::gather(MNPC1,nsd,usol[0],*eVs[0]))
      std::cerr << "*** ChorinPressCorrMixed::initElement: Detected "
		<< ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


bool ChorinPressCorr::evalInt(LocalIntegral*& elmInt, 
			      const TimeDomain& time, double detJW,
			      const Vector& N, const Matrix& dNdX,
			      const Vec3& X) const
{
  int i, j, k;
  double divU;
  double laplace, rhoDtInv;
  double B0;

  // Time coefficient
  rhoDtInv = 1.0/time.dt;
  rhoDtInv *= rho*detJW;

  // Time integration parameter
  if (fabs(time.t - time.dt) < 1.0e-8)
    B0 = 1.0;
  else
    B0 = Beta0;

  // Divergence of velocity at integration point
  divU = 0.0;
  Vector& EV = *eVs[0];
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      divU += EV((i-1)*nsd+k)*dNdX(i,k);

  // Assembling lhs terms
  if (eM && !myMats->rhsOnly) {
    Matrix& EM = *eM;

    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dNdX(i,k)*dNdX(j,k);
	laplace *= detJW;
	
	EM(j,i) += laplace;
      }
  }

  // Assembling rhs terms
  if (eS) {
    Vector& ES = *eS;

    for (i = 1;i <= N.size();i++)
      ES(i) -= B0*rhoDtInv*divU*N(i);
  }

  return getIntegralResult(elmInt);
}


bool ChorinPressCorr::evalInt(LocalIntegral*& elmInt, 
			      const TimeDomain& time,
			      double detJW,
			      const Vector& N1, 
			      const Vector& N2,
			      const Matrix& dN1dX, 
			      const Matrix& dN2dX,
			      const Vec3& X) const
{
  int i, j, k;
  double divU;
  double laplace, rhoDtInv;
  double B0;

  // Time coefficient
  rhoDtInv = 1.0/time.dt;
  rhoDtInv *= rho*detJW;

  // Time integration parameter
  if (fabs(time.t - time.dt) < 1.0e-10)
    B0 = 1.0;
  else
    B0 = Beta0;

  // Divergence of velocity at integration point
  divU = 0.0;
  Vector& EV = *eVs[0];
  for (i = 1;i <= N1.size();i++)
    for(k = 1;k <= nsd;k++)
      divU += EV((i-1)*nsd+k)*dN1dX(i,k);

  // Assembling matrix terms
  if (eM && !myMats->rhsOnly) {
    Matrix& EM = *eM;

    for (i = 1;i <= N2.size();i++)
      for (j = 1;j <= N2.size();j++) {
	laplace = 0.0;
	for (k = 1;k <= nsd;k++)
	  laplace += dN2dX(i,k)*dN2dX(j,k);
	laplace *= detJW;
	
	EM(j,i) += laplace;
      }
  }

  // Assembling rhs terms
  if (eS) {
    Vector& ES = *eS;
    
    for (i = 1;i <= N2.size();i++)
      ES(i) -= B0*rhoDtInv*divU*N2(i);
  }

  return getIntegralResult(elmInt);
}


bool ChorinPressCorr::evalBou(LocalIntegral*& elmInt, 
			      const TimeDomain& time, double detJW,
			      const Vector& N, const Matrix& dNdX,
			      const Vec3& X, const Vec3& normal) const
{
  return getIntegralResult(elmInt);
}


bool ChorinPressCorr::evalSol(Vector& s,
			      const Vector& N, const Matrix& dNdX,
			      const Vec3& X, const std::vector<int>& MNPC) const
{
  int i, k;
  Vector P;

  if (primsol.empty()) {
    std::cerr << " *** ChorinPressCorr::evalSol: No primary solution."
	      << std::endl;
    return false;
  }

  // Element pressure vector
  int ierr = 0;
  if (ierr = utl::gather(MNPC,1,primsol[0],P)) {
    std::cerr << " *** ChorinPressCorr::evalSol: Detected "
	      << ierr << " node numbers out of range." << std::endl;
    return false;
  }

  // Evaluate flux
  s.resize(nsd,true);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      s(k) -= P(i)*dNdX(i,k);

  return true;
}


bool ChorinPressCorr::evalSolScal(Vector& s,
				  const VecFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


bool ChorinPressCorr::getIntegralResult (LocalIntegral*& elmInt) const
{
  elmInt = myMats;
  return elmInt ? true : false;
}
