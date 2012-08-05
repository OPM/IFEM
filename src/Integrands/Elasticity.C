// $Id$
//==============================================================================
//!
//! \file Elasticity.C
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear and nonlinear elasticity problems.
//!
//==============================================================================

#include "Elasticity.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"
#include <iomanip>

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif


Elasticity::Elasticity (unsigned short int n, bool ax) : nsd(n), axiSymmetry(ax)
{
  if (axiSymmetry) nsd = 2;

  nDF = axiSymmetry ? 3 : nsd;
  npv = nsd; // Number of primary unknowns per node

  // Default is zero gravity
  grav[0] = grav[1] = grav[2] = 0.0;

  material = 0;
  locSys = 0;
  tracFld = 0;
  fluxFld = 0;
  bodyFld = 0;
  eM = eKm = eKg = 0;
  eS = iS = 0;
}


Elasticity::~Elasticity ()
{
  if (locSys) delete locSys;
}


void Elasticity::print (std::ostream& os) const
{
  if (axiSymmetry)
    std::cout <<"Axial-symmetric Elasticity problem\n";
  std::cout <<"Elasticity: "<< nsd <<"D, gravity =";
  for (unsigned short int d = 0; d < nsd; d++)
    std::cout <<" "<< grav[d];
  std::cout << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<Elasticity*>(this)->material = &defaultMat;
  }
  material->print(os);
}


void Elasticity::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  eM = eKm = eKg = 0;
  eS = iS = 0;

  if (mode == SIM::BUCKLING || mode == SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.clear();

  switch (mode)
    {
    case SIM::STATIC:
      eKm = 1;
      eS  = 1;
      break;

    case SIM::DYNAMIC:
      eKm = 1;
      eM  = 2;
      eS  = 1;
      break;

    case SIM::VIBRATION:
      eKm = 1;
      eM  = 2;
      break;

    case SIM::BUCKLING:
      eKm = 1;
      eKg = 2;
      break;

    case SIM::STIFF_ONLY:
      eKm = 1;
      break;

    case SIM::MASS_ONLY:
      eM = 1;
      break;

    case SIM::RHS_ONLY:
      eS = 1;
      break;

    case SIM::RECOVERY:
      maxVal.clear();
      break;

    default:
      ;
    }
}


LocalIntegral* Elasticity::getLocalIntegral (size_t nen, size_t,
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
    case SIM::BUCKLING:
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
    result->A[i].resize(nsd*nen,nsd*nen);

  if (result->b.size())
    result->b.front().resize(nsd*nen);

  return result;
}


Vec3 Elasticity::getTraction (const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


Vec3 Elasticity::getBodyforce (const Vec3& X) const
{
  Vec3 f(grav[0],grav[1],grav[2]);
  f *= material->getMassDensity(X);

  if (bodyFld)
    f += (*bodyFld)(X);

  return f;
}


bool Elasticity::haveLoads () const
{
  if (tracFld) return true;
  if (fluxFld) return true;
  if (bodyFld) return true;

  for (unsigned short int i = 0; i < nsd; i++)
    if (grav[i] != 0.0)
      if (material)
	return material->getMassDensity(Vec3()) != 0.0;

  return false;
}


void Elasticity::initIntegration (size_t, size_t nBp)
{
  tracVal.clear();
  tracVal.resize(nBp,std::make_pair(Vec3(),Vec3()));
}


bool Elasticity::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
{
  if (tracVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write boundary tractions as discrete point vectors to the VTF-file
  return vtf->writeVectors(tracVal,++nBlock,"Tractions",iStep);
}


/*!
  The strain-displacement matrix for a continuum element is formally defined as:
  \f[ \mbox{In 3D,~~}
  [B] = \left[\begin{array}{ccc}
  \frac{\partial}{\partial x} & 0 & 0 \\
  0 & \frac{\partial}{\partial y} & 0 \\
  0 & 0 & \frac{\partial}{\partial z} \\
  \frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
  0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x}
  \end{array}\right] [N] \hskip 5mm\mbox{and in 2D,~~}
  [B] = \left[\begin{array}{ccc}
  \frac{\partial}{\partial x} & 0 \\
  0 & \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [nsd][nsd*NENOD] matrix.
*/

bool Elasticity::formBmatrix (Matrix& Bmat, const Matrix& dNdX) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid dimension on dNdX, "
	      << dNdX.rows() <<"x"<< dNdX.cols() <<"."<< std::endl;
    return false;
  }

#define INDEX(i,j) i+nstrc*(j-1)

  switch (nsd) {
  case 1:

    // Strain-displacement matrix for 1D elements:
    //
    //   [B] = | d/dx | * [N]

    for (size_t i = 1; i <= nenod; i++)
      Bmat(1,i) = dNdX(i,1);
    break;

  case 2:

    // Strain-displacement matrix for 2D elements:
    //
    //         | d/dx   0   |
    //   [B] = |  0    d/dy | * [N]
    //         | d/dy  d/dx |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      // Shear strain part
      Bmat(INDEX(3,1),i) = dNdX(i,2);
      Bmat(INDEX(3,2),i) = dNdX(i,1);
    }
    break;

  case 3:

    // Strain-displacement matrix for 3D elements:
    //
    //         | d/dx   0     0   |
    //         |  0    d/dy   0   |
    //   [B] = |  0     0    d/dz | * [N]
    //         | d/dy  d/dx   0   |
    //         | d/dz   0    d/dx |
    //         |  0    d/dz  d/dy |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      Bmat(INDEX(3,3),i) = dNdX(i,3);
      // Shear strain part
      Bmat(INDEX(4,1),i) = dNdX(i,2);
      Bmat(INDEX(4,2),i) = dNdX(i,1);
      Bmat(INDEX(5,2),i) = dNdX(i,3);
      Bmat(INDEX(5,3),i) = dNdX(i,2);
      Bmat(INDEX(6,1),i) = dNdX(i,3);
      Bmat(INDEX(6,3),i) = dNdX(i,1);
    }
    break;

  default:
    std::cerr <<" *** Elasticity::formBmatrix: nsd="<< nsd << std::endl;
    return false;
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);
  return true;
}


/*!
  The strain-displacement matrix for an axially symmetric 3D continuum element
  is formally defined as:
  \f[
  [B] = \left[\begin{array}{cc}
  \frac{\partial}{\partial r} &                0            \\
                 0            & \frac{\partial}{\partial z} \\
         \frac{1}{r}          &                0            \\
  \frac{\partial}{\partial z} & \frac{\partial}{\partial r}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [2][2*NENOD] matrix.
*/

bool Elasticity::formBmatrix (Matrix& Bmat, const Vector& N, const Matrix& dNdX,
			      const double r) const
{
  const size_t nenod = N.size();
  Bmat.resize(8,nenod,true);
  if (dNdX.cols() < 2)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid dimension on dNdX, "
	      << dNdX.rows() <<"x"<< dNdX.cols() <<"."<< std::endl;
    return false;
  }
  else if (r < -epsR)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid point r < 0, "
	      << r << std::endl;
    return false;
  }

#define INDEX(i,j) i+4*(j-1)

  // Strain-displacement matrix for 3D axisymmetric elements:
  //
  //         | d/dr   0   |
  //   [B] = |  0    d/dz | * [N]
  //         | 1/r    0   |
  //         | d/dz  d/dr |

  for (size_t i = 1; i <= nenod; i++)
  {
    // Normal strain part
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    // Hoop strain part
    Bmat(INDEX(3,1),i) = r <= epsR ? dNdX(i,1) : N(i)/r;
    // Shear strain part
    Bmat(INDEX(4,1),i) = dNdX(i,2);
    Bmat(INDEX(4,2),i) = dNdX(i,1);
  }

#undef INDEX

  Bmat.resize(4,2*nenod);
  return true;
}


bool Elasticity::kinematics (const Vector& eV,
			     const Vector& N, const Matrix& dNdX, double r,
			     Matrix& B, Tensor&, SymmTensor& eps) const
{
  // Evaluate the strain-displacement matrix, B
  if (axiSymmetry)
  {
    if (!this->formBmatrix(B,N,dNdX,r))
      return false;
  }
  else
  {
    if (!this->formBmatrix(B,dNdX))
      return false;
  }

  if (eV.empty() || eps.dim() == 0)
    return true;

  // Evaluate the strains
  return B.multiply(eV,eps); // eps = B*eV
}


void Elasticity::formKG (Matrix& EM, const Vector& N, const Matrix& dNdX,
			 double r, const Tensor& sigma, double detJW) const
{
#if SP_DEBUG > 3
  std::cout <<"Elasticity::sigma =\n"<< sigma;
  std::cout <<"Elasticity::kg =";
#endif

  unsigned short int i, j;
  double kgrr = axiSymmetry && r > 0.0 ? sigma(3,3)/(r*r) : 0.0;
  for (size_t a = 1; a <= dNdX.rows(); a++)
    for (size_t b = 1; b <= dNdX.rows(); b++)
    {
      double kg = 0.0;
      for (i = 1; i <= nsd; i++)
	for (j = 1; j <= nsd; j++)
	  kg += dNdX(a,i)*sigma(i,j)*dNdX(b,j);
#if SP_DEBUG > 3
      std::cout << (b == 1 ? '\n' : ' ') << kg;
#endif

      for (i = 1; i <= nsd; i++)
	EM(nsd*(a-1)+i,nsd*(b-1)+i) += kg*detJW;

      if (kgrr > 0.0)
	EM(nsd*(a-1)+1,nsd*(b-1)+1) += N(a)*kgrr*N(b)*detJW;
    }
#if SP_DEBUG > 3
  std::cout << std::endl;
#endif
}


void Elasticity::formMassMatrix (Matrix& EM, const Vector& N,
				 const Vec3& X, double detJW) const
{
  double rhow = material->getMassDensity(X)*detJW;
  if (rhow == 0.0) return;

  for (size_t a = 1; a <= N.size(); a++)
    for (size_t b = 1; b <= N.size(); b++)
      for (unsigned short int i = 1; i <= nsd; i++)
	EM(nsd*(a-1)+i,nsd*(b-1)+i) += rhow*N(a)*N(b);
}


void Elasticity::formBodyForce (Vector& ES, const Vector& N,
				const Vec3& X, double detJW) const
{
  Vec3 f = this->getBodyforce(X);
  if (f.isZero()) return;

  f *= detJW;
  for (size_t a = 1; a <= N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += f[i-1]*N(a);
}


Vec3 Elasticity::evalSol (const Vector& eV, const Vector& N) const
{
  Vec3 u;
  if (eV.size() == N.size()*nsd)
    for (unsigned short int i = 0; i < nsd; i++)
      u[i] = eV.dot(N,i,nsd);

  return u;
}


bool Elasticity::formCinverse (Matrix& Cinv, const Vec3& X) const
{
  SymmTensor dummy(nsd,axiSymmetry); double U;
  return material->evaluate(Cinv,dummy,U,0,X,dummy,dummy,-1);
}


bool Elasticity::evalSol (Vector& s,
			  const Vector& N, const Matrix& dNdX, const Vec3& X,
			  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),eV)))
    {
      std::cerr <<" *** Elasticity::evalSol: Detected "<< ierr
		<<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the stress tensor
  if (!this->evalSol(s,eV,N,dNdX,X,true))
    return false;

  // Additional result variables?
  for (int i = 1; i <= material->getNoIntVariables(); i++)
    s.push_back(material->getInternalVariable(i));

  // Find the maximum values for each quantity
  if (maxVal.empty())
    for (size_t j = 0; j < s.size(); j++)
      maxVal.push_back(std::make_pair(X,s[j]));
  else
    for (size_t j = 0; j < s.size(); j++)
      if (fabs(s[j]) > fabs(maxVal[j].second))
	maxVal[j] = std::make_pair(X,s[j]);

  return true;
}


bool Elasticity::evalSol (Vector& s, const Vector& eV, const Vector& N,
			  const Matrix& dNdX, const Vec3& X, bool toLocal) const
{
  if (eV.size() != dNdX.rows()*nsd)
  {
    std::cerr <<" *** Elasticity::evalSol: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV.size() <<"   size(dNdX) = "
	      << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the deformation gradient, dUdX, and/or the strain tensor, eps
  Matrix Bmat;
  Tensor dUdX(nDF);
  SymmTensor eps(nsd,axiSymmetry);
  if (!this->kinematics(eV,N,dNdX,X.x,Bmat,dUdX,eps))
    return false;

  // Calculate the stress tensor through the constitutive relation
  Matrix Cmat;
  SymmTensor sigma(nsd, axiSymmetry || material->isPlaneStrain()); double U;
  if (!material->evaluate(Cmat,sigma,U,0,X,dUdX,eps))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (toLocal && locSys) sigma.transform(locSys->getTmat(X));

  s = sigma;
  if (toLocal) s.push_back(sigma.vonMises());

  return true;
}


bool Elasticity::evalSol (Vector& s, const STensorFunc& asol,
			  const Vec3& X) const
{
  s = asol(X);
  s.push_back(SymmTensor(s).vonMises());
  return true;
}


size_t Elasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd;
  else if (nsd == 2 && (axiSymmetry || material->isPlaneStrain()))
    return 5 + material->getNoIntVariables();
  else
    return nsd*(nsd+1)/2 + 1 + material->getNoIntVariables();
}


const char* Elasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i > nsd) i = 4;

  static const char* s[5] = { "u_x", "u_y", "u_z", "u_r", "displacement" };
  if (!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[axiSymmetry ? 3-i : i];

  return name.c_str();
}


const char* Elasticity::getField2Name (size_t i, const char* prefix) const
{
  size_t nStress = this->getNoFields(2);
  if (i >= nStress) return 0;

  static const char* r[4] = { "s_rr", "s_zz", "s_tt", "s_zr" };
  static const char* s[6] = { "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz" };

  static std::string name;
  if (prefix)
    name = std::string(prefix) + " ";
  else
    name.clear();

  // Number of components in the stress vector of this problem
  nStress -= 1 + material->getNoIntVariables();

  if (nsd == 1)
    name += "Axial stress";
  else if (i == 2 && nStress == 3)
    name += s[3]; // No s_zz when plane stress
  else if (i < nStress)
    name += axiSymmetry ? r[i] : s[i];
  else if (i == nStress)
    name += "von Mises stress";
  else
  {
    static char varName[32];
    material->getInternalVariable(i-nStress,varName);
    name += varName;
  }

  return name.c_str();
}


void Elasticity::printMaxVals (std::ostream& os,
                               std::streamsize precision, size_t comp) const
{
  size_t i1 = 1, i2 = maxVal.size();
  if (comp > i2)
    return;
  else if (comp > 0)
    i1 = i2 = comp;

  for (size_t i = i1; i <= i2; i++)
  {
    const char* name = this->getField2Name(i-1);
    os <<"  Max "<< name <<":";
    for (size_t j = strlen(name); j < 16; j++) std::cout <<' ';
    std::streamsize flWidth = 8 + precision;
    std::streamsize oldPrec = os.precision(precision);
    std::ios::fmtflags oldF = os.flags(std::ios::scientific | std::ios::right);
    os << std::setw(flWidth) << maxVal[i-1].second;
    os.precision(oldPrec);
    os.flags(oldF);
    std::cout <<"  X = "<< maxVal[i-1].first << std::endl;
  }
}


NormBase* Elasticity::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new ElasticityNorm(*const_cast<Elasticity*>(this),
			      asol->getStressSol());
  else
    return new ElasticityNorm(*const_cast<Elasticity*>(this));
}


ElasticityNorm::ElasticityNorm (Elasticity& p, STensorFunc* a)
  : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


size_t ElasticityNorm::getNoFields (int fld) const
{
  if (fld == 0) {
    size_t nf=1;
    for (size_t i = 0; i < prjsol.size(); i++)
      if (!prjsol.empty())
        nf++;
    return nf;
  }
  if (fld == 1)
    return anasol ? 4 : 2;

  return anasol ? 6 : 4;
}


bool ElasticityNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
			      const Vec3& X) const
{
  Elasticity& problem = static_cast<Elasticity&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCinverse(Cinv,X)) return false;

  // Evaluate the finite element stress field
  Vector sigmah, sigma, error;
  if (!problem.evalSol(sigmah,pnorm.vec.front(),fe.N,fe.dNdX,X))
    return false;

  bool planeStrain = sigmah.size() == 4 && Cinv.rows() == 3;
  if (planeStrain) sigmah.erase(sigmah.begin()+2); // Remove the sigma_zz

  double detJW = fe.detJxW;
  if (problem.isAxiSymmetric())
    detJW *= 2.0*M_PI*X.x;

  size_t ip = 0;
  // Integrate the energy norm a(u^h,u^h)
  pnorm[ip++] += sigmah.dot(Cinv*sigmah)*detJW;

  if (problem.haveLoads())
  {
    // Evaluate the body load
    Vec3 f = problem.getBodyforce(X);
    // Evaluate the displacement field
    Vec3 u = problem.evalSol(pnorm.vec.front(),fe.N);
    // Integrate the external energy (f,u^h)
    pnorm[ip] += f*u*detJW;
  }
  ip++;

  if (anasol)
  {
    // Evaluate the analytical stress field
    sigma = (*anasol)(X);
    if (sigma.size() == 4 && Cinv.rows() == 3)
      sigma.erase(sigma.begin()+2); // Remove the sigma_zz if plane strain

    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Cinv*sigma)*detJW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Cinv*error)*detJW;
  }

  size_t i, j, k;
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected stress field
      Vector sigmar(sigmah.size());
      for (j = k = 0; j < nrcmp && k < sigmar.size(); j++)
	if (!planeStrain || j != 2)
	  sigmar[k++] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Cinv*sigmar)*detJW;
      // Integrate the error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Cinv*error)*detJW;

      double l2u = sigmar.norm2();
      double l2e = error.norm2();

      // Integrate the L2-norm (sigma^r,sigma^r)
      pnorm[ip++] += l2u*l2u*detJW;
      // Integrate the error in L2-norm (sigma^r-sigma^h,sigma^r-sigma^h)
      pnorm[ip++] += l2e*l2e*detJW;

      if (anasol)
      {
	// Integrate the error in the projected solution a(u-u^r,u-u^r)
	error = sigma - sigmar;
	pnorm[ip++] += error.dot(Cinv*error)*detJW;
	ip++; // Make room for the local effectivity index here
      }
    }

  return true;
}


bool ElasticityNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			      const Vec3& X, const Vec3& normal) const
{
  Elasticity& problem = static_cast<Elasticity&>(myProblem);
  if (!problem.haveLoads()) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface traction
  Vec3 T = problem.getTraction(X,normal);
  // Evaluate the displacement field
  Vec3 u = problem.evalSol(pnorm.vec.front(),fe.N);

  double detJW = fe.detJxW;
  if (problem.isAxiSymmetric())
    detJW *= 2.0*M_PI*X.x;

  // Integrate the external energy
  pnorm[1] += T*u*detJW;
  return true;
}


bool ElasticityNorm::finalizeElement (LocalIntegral& elmInt,
				      const TimeDomain&, size_t)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 9; ip < pnorm.size(); ip += 6)
    pnorm[ip] = sqrt(pnorm[ip-4] / pnorm[3]);

  return true;
}


const char* ElasticityNorm::getName(size_t i, size_t j, const char* prefix)
{
  static const char* u[4] = {
    "a(u^h,u^h)^0.5",
    "(f,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h"
  };

  static const char* p[6] = {
    "a(u^r,u^r)^0.5",
    "a(e',e')^0.5, e'=u^r-u^h",
    "(u^r,u^r)^0.5",
    "(e',e')^0.5, e'=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index"
  };

  const char** s = j>1?p:u;

  if (!prefix)
    return s[j-1];

  static std::string name;
  name = prefix + std::string(" ");
  name += s[j-1];

  return name.c_str();
}


void ElasticityNorm::addBoundaryTerms(Vectors& gNorm, double extEnergy)
{
  gNorm[0](2) += extEnergy;
}
