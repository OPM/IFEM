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
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"


Poisson::Poisson (unsigned short int n) : nsd(n)
{
  kappa = 1.0;

  tracFld = 0;
  heatSrc = 0;
  eM = 0;
  eS = eV = 0;

  // Only the current solution is needed
  primsol.resize(1);
}


void Poisson::setMode (SIM::SolutionMode mode)
{
  myMats.rhsOnly = false;
  eM = 0;
  eS = eV = 0;

  switch (mode)
    {
    case SIM::STATIC:
      myMats.A.resize(1);
      myMats.b.resize(1);
      eM = &myMats.A[0];
      eS = &myMats.b[0];
      tracVal.clear();
      break;

    case SIM::VIBRATION:
      myMats.A.resize(1);
      eM = &myMats.A[0];
      break;

    case SIM::RHS_ONLY:
      myMats.rhsOnly = true;
      if (myMats.b.empty())
	myMats.b.resize(1);
      eS = &myMats.b[0];
      tracVal.clear();
      break;

    case SIM::RECOVERY:
      myMats.rhsOnly = true;
      if (myMats.b.empty())
	myMats.b.resize(1);
      eV = &myMats.b[0];
      break;

    default:
      myMats.A.clear();
      myMats.b.clear();
      tracVal.clear();
    }
}

double Poisson::getTraction (const Vec3& X, const Vec3& n) const
{
  if (tracFld)
    return (*tracFld)(X)*(n);
  else
    return 0;
}



bool Poisson::initElement (const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eM) eM->resize(nen,nen,true);
  if (eS) eS->resize(nen,true);

  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,1,primsol.front(),*eV)))
      std::cerr <<" *** Poisson::initElement: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;

  myMats.withLHS = true;
  return ierr == 0;
}


bool Poisson::initElementBou (const std::vector<int>& MNPC)
{
  if (eS) eS->resize(MNPC.size(),true);

  myMats.withLHS = false;
  return true;
}


bool Poisson::evalInt (LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X) const
{
  elmInt = &myMats;

  if (eM)
  {
    // Evaluate the constitutive matrix at this point
    if (!this->formCmatrix(C,X)) return false;

    CB.multiply(C,fe.dNdX,false,true).multiply(fe.detJxW); // = C*dNdX^T*|J|*w
    eM->multiply(fe.dNdX,CB,false,false,true); // EK += dNdX * CB
  }

  // Integrate heat source if defined
  if (eS && heatSrc)
    eS->add(fe.N,(*heatSrc)(X)*fe.detJxW);

  return true;
}


bool Poisson::evalBou (LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const
{
  elmInt = &myMats;
  if (!tracFld)
  {
    std::cerr <<" *** Poisson::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** Poisson::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the heat flux q at point X
  Vec3 q = (*tracFld)(X);

  // Evaluate the Neumann value -q*n
  double trac = -(q*normal);

  // Store traction value for visualization
  if (abs(trac) > 1.0e8)
    tracVal.insert(std::make_pair(X,trac*normal));

  // Integrate the Neumann value
  eS->add(fe.N,trac*fe.detJxW);

  return true;
}


bool Poisson::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
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
    std::cerr <<" *** Poisson::evalSol: No primary solution."
	      << std::endl;
    return false;
  }
  else if (!this->formCmatrix(C,X))
    return false;

  Vector Dtmp;
  int ierr = utl::gather(MNPC,1,primsol.front(),Dtmp);
  if (ierr > 0)
  {
    std::cerr <<" *** Poisson::evalSol: Detected "
	      << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

  // Evaluate the heat flux vector
  CB.multiply(C,dNdX,false,true).multiply(Dtmp,q); // q = C*dNdX^T*Dtmp
  q *= -1.0;

  return true;
}


bool Poisson::evalSol (Vector& q, const Matrix& dNdX, const Vec3& X) const
{
  if (!eV || eV->empty())
  {
    std::cerr <<" *** Poisson::evalSol: No solution vector."
	      << std::endl;
    return false;
  }
  else if (eV->size() != dNdX.rows())
  {
    std::cerr <<" *** Poisson::evalSol: Invalid solution vector."
	      <<"\n     size(eV) = "<< eV->size() <<"   size(dNdX) = "
	      << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }
  else if (!this->formCmatrix(C,X))
    return false;

  // Evaluate the heat flux vector
  CB.multiply(C,dNdX,false,true).multiply(*eV,q); // q = C*dNdX^T*eV
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


bool PoissonNorm::initElement (const std::vector<int>& MNPC)
{
  return problem.initElement(MNPC);
}


bool PoissonNorm::evalInt (LocalIntegral*& elmInt, const FiniteElement& fe,
			   const Vec3& X) const
{
  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (!eNorm)
  {
    std::cerr <<" *** PoissonNorm::evalInt: No norm array."<< std::endl;
    return false;
  }

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCmatrix(Cinv,X,true)) return false;

  // Evaluate the finite element heat flux field
  Vector sigmah;
  if (!problem.evalSol(sigmah,fe.dNdX,X)) return false;

  // Integrate the energy norm a(u^h,u^h)
  ElmNorm& pnorm = *eNorm;
  pnorm[0] += sigmah.dot(Cinv*sigmah)*fe.detJxW;
  if (anasol)
  {
    // Evaluate the analytical heat flux
    Vector sigma((*anasol)(X).ptr(),sigmah.size());
    // Integrate the energy norm a(u,u)
    pnorm[1] += sigma.dot(Cinv*sigma)*fe.detJxW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    sigma -= sigmah;
    pnorm[2] += sigma.dot(Cinv*sigma)*fe.detJxW;
  }

  return true;
}
