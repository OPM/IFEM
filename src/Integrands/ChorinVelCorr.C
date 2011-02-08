//==============================================================================
//!
//! \file ChorinVelCorr.C
//!
//! \date Sep 30 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for velocity correction in Chorin's method
//!
//==============================================================================

#include "ChorinVelCorr.h"
#include "Utilities.h"
#include "Vec3Oper.h"

ChorinVelCorr::ChorinVelCorr(unsigned short int n, double B0, 
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

  // Velocity, temporary velocity and pressure correction is needed
  int nUsols = 1;
  primsol.resize(nUsols);
  eVs.resize(nUsols);
  for (int i = 0;i < nUsols;i++)
    eVs[i] = new Vector();

  int nPsols = 1;
  psol.resize(nPsols);
  ePs.resize(nPsols);
  for (int i = 0;i < nPsols;i++)
    ePs[i] = new Vector();
}


ChorinVelCorr::~ChorinVelCorr()
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


bool ChorinVelCorr::initElement(const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eM) eM->resize(nsd*nen,nsd*nen,true);
  if (eS) eS->resize(nsd*nen,true);

  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    for (int i = 0;i < eVs.size();i++)
      if (ierr = utl::gather(MNPC,nsd,primsol[i],*eVs[i]))
	std::cerr << "*** ChorinVelCorr::initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;
  if (!ePs.empty() && !psol.empty())
    for (int i = 0;i < ePs.size();i++)
      if (ierr = utl::gather(MNPC,1,psol[i],*ePs[i]))
	std::cerr << "*** ChorinVelCorr::initElement: Detected "
		  << ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinVelCorr::initElement(const std::vector<int>& MNPC1,
				const std::vector<int>& MNPC2, 
				size_t n1)
{
  // Only velocity degrees of freedom
  const size_t nen1 = MNPC1.size();

  if (eM) eM->resize(nsd*nen1,nsd*nen1,true);
  if (eS) eS->resize(nsd*nen1,true);
  
  // Extract element velocity vector
  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    for (int i = 0;i < eVs.size();i++)
      if (ierr = utl::gather(MNPC1,nsd,primsol[i],*eVs[i]))
        std::cerr << "*** ChorinVelCorrMixed::initElement: Detected "
                  << ierr << " node numbers out of range." << std::endl;

  if (!ePs.empty() && !psol.empty())
    for (int i = 0;i < ePs.size();i++)
      if (ierr = utl::gather(MNPC2,1,psol[i],*ePs[i]))
        std::cerr << "*** ChorinVelCorrMixed::initElement: Detected "
                  << ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool ChorinVelCorr::initElementBou(const std::vector<int>& MNPC)
{				     
  const size_t nen = MNPC.size();
  
  if (eS) eS->resize(nsd*nen,true);
  
  myMats->withLHS = false;
  return true;
}


bool ChorinVelCorr::initElementBou(const std::vector<int>& MNPC1,
				   const std::vector<int>& MNPC2, 
				   size_t n1)
{
  // Only velocity degrees of freedom
  const size_t nen1 = MNPC1.size();

  if (eS) eS->resize(nsd*nen1,true);

  // Extract element velocity vector
  int ierr = 0;
  if (!eVs.empty() && !primsol.empty())
    if (ierr = utl::gather(MNPC1,nsd,primsol[0],*eVs[0]))
      std::cerr << "*** ChorinVelCorrMixed:initElementBou: Detected "
                << ierr << " node numbers out of range." << std::endl;

  // Extract element pressure vector
  if (!ePs.empty() && !psol.empty())
    if (ierr = utl::gather(MNPC2,1,psol[0],*ePs[0]))
      std::cerr << "*** ChorinVelCorr:initElementBou: Detected "
                << ierr << " node numbers out of range." << std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


bool ChorinVelCorr::evalInt(LocalIntegral*& elmInt, 
			    const TimeDomain& time, double detJW,
			    const Vector& N, const Matrix& dNdX,
			    const Vec3& X) const
{
  int i, j, k;
  double divU, P;
  double B0;
  double mass, rhoDtInv;
  Vector U(nsd), dPdX(nsd);

  // Time coefficient
  rhoDtInv = 1.0/time.dt;
  rhoDtInv *= rho*detJW;

  // Time integration parameter
  if (fabs(time.t-time.dt) < 1.0e-8)
    B0 = 1.0;
  else
    B0 = Beta0;

  // Velocity at integration point
  U.fill(0.0);
  Vector& EV = *eVs[0];
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      U(k) += EV((i-1)*nsd + k)*N(i);

  // Pressure at integration point
  Vector& EP = *ePs[0];
  P = 0.0;
  for (i = 1;i <= N.size();i++)
    P += EP(i)*N(i);

  // Pressure gradient at integration point
  dPdX.fill(0.0);
  for (i = 1;i <= N.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdX(k) += EP(i)*dNdX(i,k);

  // Assembling lhs terms
  if (eM) {
    Matrix& EM = *eM;

    for (i = 1;i <= N.size();i++)
      for (j = 1;j <= N.size();j++) {
	mass = N(i)*N(j)*B0*rhoDtInv;

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd+k,(i-1)*nsd+k) += mass;
      }
  }

  // Assembling rhs terms
  if (eS) {
    Vector& ES = *eS;
    
    for (i = 1;i <= N.size();i++)
      for (k = 1;k <= nsd;k++)
	//ES((i-1)*nsd+k) += B0*rhoDtInv*U(k)*N(i) - dPdX(k)*N(i)*detJW;
	ES((i-1)*nsd+k) += B0*rhoDtInv*U(k)*N(i) + P*dNdX(i,k)*detJW;
  }

  return getIntegralResult(elmInt);
}


bool ChorinVelCorr::evalInt(LocalIntegral*& elmInt,
			    const TimeDomain& time, double detJW,
			    const Vector& N1, const Vector& N2,
			    const Matrix& dN1dX, const Matrix& dN2dX,
			    const Vec3& X) const
{
  int i, j, k;
  double divU, P;
  double B0;
  double mass, rhoDtInv;
  Vector U(nsd), dPdX(nsd);

 
  // Time coefficient
  rhoDtInv = 1.0/time.dt;
  rhoDtInv *= rho*detJW;

  // Time integration parameter
  if (fabs(time.t-time.dt) < 1.0e-10)
    B0 = 1.0;
  else
    B0 = Beta0;

  // Velocity at integration point
  U.fill(0.0);
  Vector& EV = *eVs[0];
  for (i = 1;i <= N1.size();i++)
    for (k = 1;k <= nsd;k++)
      U(k) += EV((i-1)*nsd + k)*N1(i);

  // Pressure at integration point
  Vector& EP = *ePs[0];
  P = 0.0;
  for (i = 1;i <= N2.size();i++)
    P += EP(i)*N2(i);

  // Pressure gradient at integration point
  dPdX.fill(0.0);
  for (i = 1;i <= N2.size();i++)
    for (k = 1;k <= nsd;k++)
      dPdX(k) += EP(i)*dN2dX(i,k);

  // Assembling lhs terms
  if (eM) {
    Matrix& EM = *eM;

    for (i = 1;i <= N1.size();i++)
      for (j = 1;j <= N1.size();j++) {
	mass = N1(i)*N1(j)*B0*rhoDtInv;

	for (k = 1;k <= nsd;k++)
	  EM((j-1)*nsd+k,(i-1)*nsd+k) += mass;
      }
  }

  // Assembling rhs terms
  if (eS) {
    Vector& ES = *eS;
    
    for (i = 1;i <= N1.size();i++)
      for (k = 1;k <= nsd;k++)
	//ES((i-1)*nsd+k) += B0*rhoDtInv*U(k)*N1(i) - dPdX(k)*N1(i)*detJW;
	ES((i-1)*nsd+k) += B0*rhoDtInv*U(k)*N1(i) + P*dN1dX(i,k)*detJW;
  }

  return getIntegralResult(elmInt);  
}


bool ChorinVelCorr::evalBou(LocalIntegral*& elmInt, 
			    const TimeDomain& time, double detJW,
			    const Vector& N, const Matrix& dNdX,
			    const Vec3& X, const Vec3& normal) const
{
  if (incPressure)
    return getIntegralResult(elmInt);
  
  if (!tracFld) {
      std::cerr <<" *** ChorinVelCorr::evalBou: No tractions."<< std::endl;
      return false;
  }
  else if (!eS) {
    std::cerr <<" *** ChorinVelCorr::evalBou: No load vector."<< std::endl;
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

  return getIntegralResult(elmInt);
}


bool ChorinVelCorr::evalSol(Vector& s,
			    const Vector& N, const Matrix& dNdX,
			    const Vec3& X, const std::vector<int>& MNPC) const
{
  return false;
}


bool ChorinVelCorr::evalSolScal(Vector& s,
				  const VecFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


bool ChorinVelCorr::getIntegralResult (LocalIntegral*& elmInt) const
{
  elmInt = myMats;
  return elmInt ? true : false;
}
