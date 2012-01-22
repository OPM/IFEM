// $Id$
//==============================================================================
//!
//! \file KirchhoffLovePlate.C
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for linear Kirchhoff-Love thin plate problems.
//!
//==============================================================================

#include "KirchhoffLovePlate.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"


KirchhoffLovePlate::KirchhoffLovePlate (unsigned short int n) : nsd(n)
{
  npv = 1; // Number of primary unknowns per node

  gravity = 0.0;
  thickness = 0.1;

  material = 0;
  locSys = 0;
  presFld = 0;
  eM = eK = 0;
  eS = 0;

  // Only the current solution is needed
  primsol.resize(1);
}


KirchhoffLovePlate::~KirchhoffLovePlate ()
{
  if (locSys) delete locSys;
}


void KirchhoffLovePlate::print (std::ostream& os) const
{
  std::cout <<"KirchhoffLovePlate: thickness = "<< thickness
	    <<", gravity = "<< gravity << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLovePlate*>(this)->material = &defaultMat;
  }

  material->print(os);
}


void KirchhoffLovePlate::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  eM = eK = 0;
  eS = 0;

  switch (mode)
    {
    case SIM::STATIC:
      eK = 1;
      eS = 1;
      presVal.clear();
      break;

    case SIM::DYNAMIC:
      eK = 1;
      eM = 2;
      eS = 1;
      presVal.clear();
      break;

    case SIM::VIBRATION:
      eK = 1;
      eM = 2;
      break;

    case SIM::STIFF_ONLY:
      eK = 1;
      break;

    case SIM::MASS_ONLY:
      eM = 1;
      break;

    case SIM::RHS_ONLY:
      eS = 1;
      presVal.clear();
      break;

    case SIM::RECOVERY:
      break;

    default:
      presVal.clear();
    }
}


LocalIntegral* KirchhoffLovePlate::getLocalIntegral (size_t nen, size_t,
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

    case SIM::DYNAMIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann?0:2,1);
      break;

    case SIM::VIBRATION:
      result->withLHS = true;
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
    case SIM::MASS_ONLY:
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
      ;
    }

  for (size_t i = 0; i < result->A.size(); i++)
    result->A[i].resize(nen,nen);

  if (result->b.size())
    result->b.front().resize(nen);

  return result;
}


double KirchhoffLovePlate::getPressure (const Vec3& X) const
{
  double p = material->getMassDensity(X)*gravity*thickness;

  if (presFld)
    p += (*presFld)(X);

  return p;
}


bool KirchhoffLovePlate::haveLoads () const
{
  if (presFld) return true;

  if (gravity != 0.0 && material)
    return material->getMassDensity(Vec3()) != 0.0;

  return false;
}


bool KirchhoffLovePlate::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
{
  if (presVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write surface pressures as discrete point vectors to the VTF-file
  if (!vtf->writeVectors(presVal,++nBlock))
    return false;

  return vtf->writeVblk(nBlock,"Pressure",1,iStep);
}


/*!
  The strain-displacement matrix for a Kirchhoff-Love plate element is formally
  defined as:
  \f[
  [B] = \left[\begin{array}{c}
   \frac{\partial^2}{\partial x^2} \\
   \frac{\partial^2}{\partial y^2} \\
  2\frac{\partial^2}{\partial x\partial y}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [NENOD] row-vector.
*/

bool KirchhoffLovePlate::formBmatrix (Matrix& Bmat,
				      const Matrix3D& d2NdX2) const
{
  const size_t nenod = d2NdX2.dim(1);
  const size_t nstrc = nsd*(nsd+1)/2;

  Bmat.resize(nstrc,nenod,true);

  if (d2NdX2.dim(2) != nsd || d2NdX2.dim(3) != nsd)
  {
    std::cerr <<" *** KirchhoffLovePlate::formBmatrix: Invalid dimension on"
	      <<" d2NdX2, "<< d2NdX2.dim(1) <<"x"<< d2NdX2.dim(2)
	      <<"x"<< d2NdX2.dim(3) <<"."<< std::endl;
    return false;
  }

  for (size_t i = 1; i <= nenod; i++)
    if (nsd == 1)
      Bmat(1,i) = d2NdX2(i,1,1);
    else
    {
      Bmat(1,i) = d2NdX2(i,1,1);
      Bmat(2,i) = d2NdX2(i,2,2);
      Bmat(3,i) = d2NdX2(i,1,2)*2.0;
    }

  return true;
}


bool KirchhoffLovePlate::formCmatrix (Matrix& C, const Vec3& X,
				      bool invers) const
{
  SymmTensor dummy(nsd); double U;
  if (!material->evaluate(C,dummy,U,0,X,dummy,dummy, invers ? -1 : 1))
    return false;

  double factor = thickness*thickness*thickness/12.0;
  C.multiply(invers ? 1.0/factor : factor);
  return true;
}


void KirchhoffLovePlate::formMassMatrix (Matrix& EM, const Vector& N,
					 const Vec3& X, double detJW) const
{
  double rho = material->getMassDensity(X)*thickness;

  if (rho != 0.0)
    EM.outer_product(N,N*rho*detJW,true);
}


void KirchhoffLovePlate::formBodyForce (Vector& ES, const Vector& N,
					const Vec3& X, double detJW) const
{
  double p = this->getPressure(X);
  if (p != 0.0)
  {
    ES.add(N,p*detJW);
    presVal[X].z = p; // Store pressure value for visualization
  }
}


bool KirchhoffLovePlate::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eK)
  {
    // Compute the strain-displacement matrix B from d2NdX2
    Matrix Bmat;
    if (!this->formBmatrix(Bmat,fe.d2NdX2)) return false;

    // Evaluate the constitutive matrix at this point
    Matrix Cmat;
    if (!this->formCmatrix(Cmat,X)) return false;

    // Integrate the stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(fe.detJxW); // CB = C*B*|J|*w
    elMat.A[eK-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);

  return true;
}


bool KirchhoffLovePlate::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& normal) const
{
  std::cerr <<" *** KirchhoffLovePlate::evalBou not implemented."<< std::endl;
  return false;
}


bool KirchhoffLovePlate::evalSol (Vector& s, const Vector&, const Matrix&,
				  const Matrix3D& d2NdX2, const Vec3& X,
				  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  int ierr = 0;
  Vector eV;
  if (!primsol.front().empty())
    if ((ierr = utl::gather(MNPC,1,primsol.front(),eV)))
    {
      std::cerr <<" *** KirchhoffLovePlate::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the stress resultant tensor
  return this->evalSol(s,eV,d2NdX2,X,true);
}


bool KirchhoffLovePlate::evalSol (Vector& s, const Vector& eV,
				  const Matrix3D& d2NdX2, const Vec3& X,
				  bool toLocal) const
{
  if (eV.empty())
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: No displacement vector."
	      << std::endl;
    return false;
  }
  else if (eV.size() != d2NdX2.dim(1))
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV.size() <<"   size(d2NdX2) = "
	      << d2NdX2.dim(1) <<","<< d2NdX2.dim(2)*d2NdX2.dim(3) << std::endl;
    return false;
  }

  // Compute the strain-displacement matrix B from d2NdX2
  Matrix Bmat;
  if (!this->formBmatrix(Bmat,d2NdX2))
    return false;

  // Evaluate the constitutive matrix at this point
  Matrix Cmat;
  if (!this->formCmatrix(Cmat,X))
    return false;

  // Evaluate the curvature tensor
  SymmTensor kappa(nsd), m(nsd);
  if (!Bmat.multiply(eV,kappa)) // kappa = B*eV
    return false;

  // Evaluate the stress resultant tensor
  if (!Cmat.multiply(-1.0*kappa,m)) // m = -C*kappa
    return false;

  // Congruence transformation to local coordinate system at current point
  if (toLocal && locSys) m.transform(locSys->getTmat(X));

  s = m;
  return true;
}


bool KirchhoffLovePlate::evalSol (Vector& s, const STensorFunc& asol,
				  const Vec3& X) const
{
  s = asol(X);
  return true;
}


size_t KirchhoffLovePlate::getNoFields (int fld) const
{
  return fld < 2 ? 1 : nsd*(nsd+1)/2;
}


const char* KirchhoffLovePlate::getField1Name (size_t,
					       const char* prefix) const
{
  if (!prefix) return "w";

  static std::string name;
  name = prefix + std::string(" w");

  return name.c_str();
}


const char* KirchhoffLovePlate::getField2Name (size_t i,
					       const char* prefix) const
{
  if (i >= 3) return 0;

  static const char* s[6] = { "m_xx", "m_yy", "m_xy" };

  static std::string name;
  if (prefix)
    name = std::string(prefix) + " ";
  else
    name.clear();

  name += s[i];

  return name.c_str();
}


NormBase* KirchhoffLovePlate::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this),
				      asol->getStressSol());
  else
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this));
}


KirchhoffLovePlateNorm::KirchhoffLovePlateNorm (KirchhoffLovePlate& p,
						STensorFunc* a)
  : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


size_t KirchhoffLovePlateNorm::getNoFields () const
{
  size_t nf = anasol ? 4 : 2;
  for (size_t i = 0; i < prjsol.size(); i++)
    if (!prjsol.empty())
       nf += anasol ? 3 : 2;

  return nf;
}


bool KirchhoffLovePlateNorm::evalInt (LocalIntegral& elmInt,
				      const FiniteElement& fe,
				      const Vec3& X) const
{
  KirchhoffLovePlate& problem = static_cast<KirchhoffLovePlate&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCmatrix(Cinv,X,true)) return false;

  // Evaluate the finite element stress field
  Vector mh, m, error;
  if (!problem.evalSol(mh,pnorm.vec.front(),fe.d2NdX2,X))
    return false;

  // Integrate the energy norm a(w^h,w^h)
  pnorm[0] += mh.dot(Cinv*mh)*fe.detJxW;
  if (problem.haveLoads())
  {
    // Evaluate the body load
    double p = problem.getPressure(X);
    // Evaluate the displacement field
    double w = pnorm.vec.front().dot(fe.N);
    // Integrate the external energy (p,w^h)
    pnorm[1] += p*w*fe.detJxW;
  }

  size_t ip = 2;
  if (anasol)
  {
    // Evaluate the analytical stress resultant field
    m = (*anasol)(X);

    // Integrate the energy norm a(w,w)
    pnorm[ip++] += m.dot(Cinv*m)*fe.detJxW;
    // Integrate the error in energy norm a(w-w^h,w-w^h)
    error = m - mh;
    pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
  }

  size_t i, j;
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected stress field
      Vector mr(mh.size());
      for (j = 0; j < nrcmp; j++)
	mr[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(w^r,w^r)
      pnorm[ip++] += mr.dot(Cinv*mr)*fe.detJxW;
      // Integrate the error in energy norm a(w^r-w^h,w^r-w^h)
      error = mr - mh;
      pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

      if (anasol)
      {
        // Integrate the error in the projected solution a(w-w^r,w-w^r)
        error = m - mr;
        pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
      }
    }

  return true;
}


bool KirchhoffLovePlateNorm::evalBou (LocalIntegral& elmInt,
				      const FiniteElement& fe,
				      const Vec3& X, const Vec3& normal) const
{
  std::cerr <<" *** KirchhoffLovePlateNorm::evalBou not included."<< std::endl;
  return false;
}
