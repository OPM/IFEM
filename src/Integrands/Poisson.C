// $Id$
//==============================================================================
//!
//! \file Poisson.C
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Integrand implementations for Poisson problems.
//!
//==============================================================================

#include "Poisson.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"


Poisson::Poisson (unsigned short int n) : nsd(n)
{
  npv = 1; // One primary unknown per node (scalar equation)

  kappa = 1.0;

  tracFld = 0;
  fluxFld = 0;
  heatSrc = 0;

  // Only the current solution is needed
  primsol.resize(1);
}


double Poisson::getHeat (const Vec3& X) const
{
  return heatSrc ? (*heatSrc)(X) : 0.0;
}


double Poisson::getTraction (const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X)*n;
  else
    return 0.0;
}


LocalIntegral* Poisson::getLocalIntegral (size_t nen, size_t,
                                          bool neumann) const
{
  ElmMats* result = new ElmMats;
  switch (m_mode)
  {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:1,1);
      break;

    case SIM::VIBRATION:
      result->withLHS = true;
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->rhsOnly = true;
      result->resize(0,1);
      break;

    case SIM::RECOVERY:
      result->rhsOnly = true;
      break;

    default:
      break;
  }

  if (result->A.size())
    result->A.front().resize(nen,nen);

  if (result->b.size())
    result->b.front().resize(nen);

  return result;
}


bool Poisson::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A.empty())
  {
    // Evaluate the constitutive matrix at this point
    Matrix C, CB;
    this->formCmatrix(C,X);

    // Integrate the coefficient matrix
    CB.multiply(C,fe.dNdX,false,true).multiply(fe.detJxW); // = C*dNdX^T*|J|*w
    elMat.A.front().multiply(fe.dNdX,CB,false,false,true); // EK += dNdX * CB
  }

  // Integrate heat source, if defined
  if (heatSrc && !elMat.b.empty())
    elMat.b.front().add(fe.N,(*heatSrc)(X)*fe.detJxW);

  return true;
}


bool Poisson::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** Poisson::evalBou: No tractions."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty())
  {
    std::cerr <<" *** Poisson::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the Neumann value -q(X)*n
  double trac = -this->getTraction(X,normal);

  // Store traction value for visualization
  if (fe.iGP < tracVal.size() && abs(trac) > 1.0e-8)
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += trac*normal;
  }

  // Integrate the Neumann value
  elMat.b.front().add(fe.N,trac*fe.detJxW);

  return true;
}


void Poisson::initIntegration (size_t, size_t nBp)
{
  tracVal.clear();
  tracVal.resize(nBp,std::make_pair(Vec3(),Vec3()));
}


bool Poisson::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
{
  if (tracVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write boundary tractions as discrete point vectors to the VTF-file
  return vtf->writeVectors(tracVal,++nBlock,"Tractions",iStep);
}


bool Poisson::formCmatrix (Matrix& C, const Vec3&, bool invers) const
{
  C.diag(invers && kappa != 0.0 ? 1.0/kappa : kappa, nsd);
  return true;
}


bool Poisson::evalSol (Vector& q, const Vector&,
		       const Matrix& dNdX, const Vec3& X,
		       const std::vector<int>& MNPC) const
{
  if (primsol.front().empty())
  {
    std::cerr <<" *** Poisson::evalSol: No primary solution."<< std::endl;
    return false;
  }

  Vector eV;
  int ierr = utl::gather(MNPC,1,primsol.front(),eV);
  if (ierr > 0)
  {
    std::cerr <<" *** Poisson::evalSol: Detected "<< ierr
	      <<" node numbers out of range."<< std::endl;
    return false;
  }

  return this->evalSol(q,eV,dNdX,X);
}


bool Poisson::evalSol (Vector& q, const Vector& eV,
                       const Matrix& dNdX, const Vec3& X) const
{
  if (eV.size() != dNdX.rows())
  {
    std::cerr <<" *** Poisson::evalSol: Invalid solution vector."
              <<"\n     size(eV) = "<< eV.size() <<"   size(dNdX) = "
              << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the constitutive matrix at this point
  Matrix C, CB;
  this->formCmatrix(C,X);

  // Evaluate the heat flux vector
  CB.multiply(C,dNdX,false,true).multiply(eV,q); // q = C*dNdX^T*eV
  q *= -1.0;

  return true;
}


bool Poisson::evalSol (Vector& s, const VecFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


const char* Poisson::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "u";

  static std::string name;
  name = prefix + std::string(" u");

  return name.c_str();
}


const char* Poisson::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return 0;

  static const char* s[3] = { "q_x","q_y","q_z" };
  if (!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


NormBase* Poisson::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new PoissonNorm(*const_cast<Poisson*>(this),
			   asol->getScalarSecSol());
  else
    return new PoissonNorm(*const_cast<Poisson*>(this));
}


PoissonNorm::PoissonNorm (Poisson& p, VecFunc* a) : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


size_t PoissonNorm::getNoFields () const
{
  size_t nf = anasol ? 4 : 2;
  for (size_t i = 0; i < prjsol.size(); i++)
    if (!prjsol.empty())
      nf += anasol ? 4 : 2;

  return nf;
}


bool PoissonNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  problem.formCmatrix(Cinv,X,true);

  // Evaluate the finite element heat flux field
  Vector sigmah, sigma, error;
  if (!problem.evalSol(sigmah,pnorm.vec.front(),fe.dNdX,X))
    return false;

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += sigmah.dot(Cinv*sigmah)*fe.detJxW;
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (h,u^h)
  pnorm[1] += problem.getHeat(X)*u*fe.detJxW;

  size_t ip = 2;
  if (anasol)
  {
    // Evaluate the analytical heat flux
    sigma.fill((*anasol)(X).ptr(),nrcmp);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Cinv*sigma)*fe.detJxW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
  }

  size_t i, j;
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected heat flux field
      Vector sigmar(nrcmp);
      for (j = 0; j < nrcmp; j++)
	sigmar[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Cinv*sigmar)*fe.detJxW;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

      if (anasol)
      {
	// Integrate the error in the projected solution a(u-u^r,u-u^r)
	error = sigma - sigmar;
	pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
	ip++; // Make room for the local effectivity index here
      }
    }

  return true;
}


bool PoissonNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X, const Vec3& normal) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface heat flux
  double t = problem.getTraction(X,normal);
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);

  // Integrate the external energy (t,u^h)
  pnorm[1] += t*u*fe.detJxW;

  return true;
}


bool PoissonNorm::finalizeElement (LocalIntegral& elmInt,
				   const TimeDomain&, size_t)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 5; ip+2 < pnorm.size(); ip += 4)
    pnorm[ip+2] = sqrt(pnorm[ip] / pnorm[3]);

  return true;
}
