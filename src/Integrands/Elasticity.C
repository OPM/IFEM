// $Id: Elasticity.C,v 1.1 2011-02-08 09:06:02 kmo Exp $
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
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "VTF.h"


Elasticity::Elasticity (unsigned short int n, bool ps)
{
  nsd = n;
  planeStress = ps;

  // Default material properties - typical values for steel (SI units)
  Emod = 2.05e11;
  nu = 0.29;
  rho = 7.85e3;

  myMats = new ElmMats();

  locSys = 0;
  tracFld = 0;
  bodyFld = 0;
  eM = eKm = eKg = 0;
  eS = eV = iS = 0;
}


Elasticity::~Elasticity ()
{
  if (myMats) delete myMats;
  if (locSys) delete locSys;
}


void Elasticity::print (std::ostream& os) const
{
  std::cout <<"Elasticity: "<< nsd <<"D";
  if (nsd == 2)
    std::cout <<" plane "<< (planeStress ? "stress" : "strain");
  std::cout <<", E = "<< Emod <<", nu = "<< nu <<", rho = "<< rho
	    <<", gravity =";
  for (unsigned short int d = 0; d < nsd; d++)
    std::cout <<" "<< g[d];
  std::cout << std::endl;
}


void Elasticity::setMode (SIM::SolutionMode mode)
{
  if (!myMats) return;

  myMats->rhsOnly = false;
  eM = eKm = eKg = 0;
  eS = eV = iS = 0;

  switch (mode)
    {
    case SIM::STATIC:
      myMats->resize(1,1);
      eKm = &myMats->A[0];
      eS  = &myMats->b[0];
      tracVal.clear();
      break;

    case SIM::DYNAMIC:
      myMats->resize(2,1);
      eKm = &myMats->A[0];
      eM  = &myMats->A[1];
      eS  = &myMats->b[0];
      tracVal.clear();
      break;

    case SIM::VIBRATION:
      myMats->resize(2,0);
      eKm = &myMats->A[0];
      eM  = &myMats->A[1];
      break;

    case SIM::BUCKLING:
      myMats->resize(2,1);
      eKm = &myMats->A[0];
      eKg = &myMats->A[1];
      eV  = &myMats->b[0];
      break;

    case SIM::STIFF_ONLY:
      myMats->resize(1,0);
      eKm = &myMats->A[0];
      break;

    case SIM::MASS_ONLY:
      myMats->resize(1,0);
      eM  = &myMats->A[0];
      break;

    case SIM::RHS_ONLY:
      myMats->rhsOnly = true;
      if (myMats->b.empty())
	myMats->b.resize(1);
      eS = &myMats->b[0];
      tracVal.clear();
      break;

    case SIM::RECOVERY:
      myMats->rhsOnly = true;
      if (myMats->b.empty())
	myMats->b.resize(1);
      eV = &myMats->b[0];
      break;

    default:
      myMats->resize(0,0);
      tracVal.clear();
    }
}


bool Elasticity::getIntegralResult (LocalIntegral*& elmInt) const
{
  elmInt = myMats;
  return elmInt ? true : false;
}


void Elasticity::setTraction (TractionFunc* tf)
{
  tracFld = tf;

  if (myMats)
    myMats->rhsOnly = true;
}


bool Elasticity::initElement (const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eKm) eKm->resize(nsd*nen,nsd*nen,true);
  if (eKg) eKg->resize(nsd*nen,nsd*nen,true);
  if (eM)  eM->resize(nsd*nen,nsd*nen,true);
  if (eS)  eS->resize(nsd*nen,true);

  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),*eV)))
      std::cerr <<" *** Elasticity::initElement: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = true;
  return ierr == 0;
}


bool Elasticity::initElementBou (const std::vector<int>& MNPC)
{
  const size_t nen = MNPC.size();

  if (eS) eS->resize(nsd*nen,true);

  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),*eV)))
      std::cerr <<" *** Elasticity::initElement: Detected "
                << ierr <<" node numbers out of range."<< std::endl;

  myMats->withLHS = false;
  return ierr == 0;
}


Vec3 Elasticity::getTraction (const Vec3& X, const Vec3& n) const
{
  if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


Vec3 Elasticity::getBodyforce (const Vec3& X) const
{
  Vec3 f(g[0],g[1],g[2]);
  f *= this->getMassDensity(X);

  if (bodyFld)
    f += (*bodyFld)(X);

  return f;
}


bool Elasticity::haveLoads () const
{
  if (tracFld) return true;
  if (bodyFld) return true;

  if (rho != 0.0)
    for (unsigned short int i = 0; i < nsd; i++)
      if (g[i] != 0.0)
	return true;

  return false;
}


bool Elasticity::writeGlvT (VTF* vtf, int iStep, int& nBlock) const
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

bool Elasticity::kinematics (const Matrix& dNdX, SymmTensor& eps) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr <<" *** Elasticity::kinematics: Invalid dimension on dNdX, "
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
    std::cerr <<" *** Elasticity::kinematics: nsd="<< nsd << std::endl;
    return false;
  }

  Bmat.resize(nstrc,nsd*nenod);

  // Evaluate the strains
  if (eV && !eV->empty() && eps.dim() > 0)
    return Bmat.multiply(*eV,eps); // eps = B*eV

  return true;
}


bool Elasticity::formBmatrix (const Matrix& dNdX) const
{
  static SymmTensor dummy(0);
  return this->Elasticity::kinematics(dNdX,dummy);
}


/*!
  The consitutive matrix for Isotrophic linear elastic problems
  is defined as follows:

  For 2D plain stress: \f[
  [C] = \frac{E}{(1-\nu^2)} \left[\begin{array}{ccc}
  1 & \ \ \nu & 0 \\
  \nu & \ \ 1 & 0 \\
  0 & \ \ 0 & \frac{1}{2}(1-\nu)
  \end{array}\right] \f]

  For 2D plain strain: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{ccc}
  1-\nu & \nu & 0 \\
  \nu & 1-\nu & 0 \\
  0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For 3D: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccccc}
  1-\nu & \nu & \nu & 0 & 0 & 0 \\
  \nu & 1-\nu & \nu & 0 & 0 & 0 \\
  \nu & \nu & 1-\nu & 0 & 0 & 0 \\
  0 & 0 & 0 & \frac{1}{2}-\nu & 0 & 0 \\
  0 & 0 & 0 & 0 & \frac{1}{2}-\nu & 0 \\
  0 & 0 & 0 & 0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]
*/

bool Elasticity::formCmatrix (Matrix& C, const Vec3&, bool invers) const
{
  const double one = 1.0;
  const size_t nst = nsd*(nsd+1)/2;
  C.resize(nst,nst,true);
  if (nsd == 1)
  {
    // Special for 1D problems
    C(1,1) = invers ? one/Emod : Emod;
    return true;
  }
  else if (nu < 0.0 || nu >= 0.5)
  {
    std::cerr <<" *** Elasticity::formCmatrix: Poisson's ratio "<< nu
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  if (invers)
    if (nsd == 3 || (nsd == 2 && planeStress))
    {
      C(1,1) = one / Emod;
      C(2,1) = -nu / Emod;
    }
    else // 2D plain strain
    {
      C(1,1) = (one - nu*nu) / Emod;
      C(2,1) = (-nu - nu*nu) / Emod;
    }

  else
    if (nsd == 2 && planeStress)
    {
      C(1,1) = Emod / (one - nu*nu);
      C(2,1) = C(1,1) * nu;
    }
    else // 2D plain strain or 3D
    {
      double fact = Emod / ((one + nu) * (one - nu - nu));
      C(1,1) = fact * (one - nu);
      C(2,1) = fact * nu;
    }

  C(1,2) = C(2,1);
  C(2,2) = C(1,1);

  const double G = Emod / (2.0 + nu + nu);
  C(nsd+1,nsd+1) = invers ? one / G : G;

  if (nsd > 2)
  {
    C(3,1) = C(2,1);
    C(3,2) = C(2,1);
    C(1,3) = C(2,1);
    C(2,3) = C(2,1);
    C(3,3) = C(1,1);
    C(5,5) = C(4,4);
    C(6,6) = C(4,4);
  }

  return true;
}


bool Elasticity::constitutive (Matrix& C, SymmTensor& sigma, double&,
			       const SymmTensor& eps, const Vec3& X,
			       char calcStress) const
{
  // Set up the constitutive matrix, C, at this point
  if (!this->formCmatrix(C,X))
    return false;
  else if (!calcStress)
    return true;

  // Calculate the stress vector, sigma = C*eps
  return C.multiply(eps,sigma);
}


void Elasticity::formKG (Matrix& EM, const Matrix& dNdX,
			 const Tensor& sigma, double detJW) const
{
#ifdef INT_DEBUG
  std::cout <<"Elasticity::eV ="<< *eV;
  std::cout <<"Elasticity::sigma =\n"<< sigma;
  std::cout <<"Elasticity::kg =";
#endif

  unsigned short int i, j;
  for (size_t a = 1; a <= dNdX.rows(); a++)
    for (size_t b = 1; b <= dNdX.rows(); b++)
    {
      double kg = 0.0;
      for (i = 1; i <= nsd; i++)
	for (j = 1; j <= nsd; j++)
	  kg += dNdX(a,i)*sigma(i,j)*dNdX(b,j);
#ifdef INT_DEBUG
      std::cout << (b == 1 ? '\n' : ' ') << kg;
#endif

      for (i = 1; i <= nsd; i++)
	EM(nsd*(a-1)+i,nsd*(b-1)+i) += kg*detJW;
    }
#ifdef INT_DEBUG
  std::cout << std::endl;
#endif
}


void Elasticity::formMassMatrix (Matrix& EM, const Vector& N,
				 const Vec3& X, double detJW) const
{
  double rhow = this->getMassDensity(X)*detJW;
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


Vec3 Elasticity::evalSol (const Vector& N) const
{
  Vec3 u;
  if (eV && eV->size() == N.size()*nsd)
    for (unsigned short int i = 0; i < nsd; i++)
      u[i] = eV->dot(N,i,nsd);

  return u;
}


bool Elasticity::evalSol (Vector& s, const Vector&,
			  const Matrix& dNdX, const Vec3& X,
			  const std::vector<int>& MNPC) const
{
  int ierr = 0;
  if (eV && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,nsd,primsol.front(),*eV)))
    {
      std::cerr <<" *** Elasticity::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the strain tensor, eps
  SymmTensor eps(nsd);
  if (!this->kinematics(dNdX,eps))
    return false;

  // Calculate the stress tensor through the constitutive relation
  SymmTensor sigma(nsd); double U = 0.0;
  if (!this->constitutive(Cmat,sigma,U,eps,X))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (locSys) sigma.transform(locSys->getTmat(X));

  s = sigma;
  s.push_back(sigma.vonMises());
  return true;
}


bool Elasticity::evalSol (Vector& s, const Matrix& dNdX, const Vec3& X) const
{
  if (!eV || eV->empty())
  {
    std::cerr <<" *** Elasticity::evalSol: No displacement vector."
	      << std::endl;
    return false;
  }
  else if (eV->size() != dNdX.rows()*nsd)
  {
    std::cerr <<" *** Elasticity::evalSol: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV->size() <<"   size(dNdX) = "
	      << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the strain tensor, eps
  SymmTensor eps(nsd);
  if (!this->kinematics(dNdX,eps))
    return false;

  // Calculate the stress tensor through the constitutive relation
  SymmTensor sigma(nsd); double U = 0.0;
  if (!this->constitutive(Cmat,sigma,U,eps,X))
    return false;

  s = sigma;
  return true;
}


bool Elasticity::evalSol (Vector& s, const TensorFunc& sol, const Vec3& X) const
{
  s = sol(X);
  s.push_back(SymmTensor(s).vonMises());
  return true;
}


size_t Elasticity::getNoFields () const
{
  return nsd*(nsd+1)/2 + 1;
}


const char* Elasticity::getFieldLabel (size_t i, const char* prefix) const
{
  static const char* s1[1] = { "Axial stress" };
  static const char* s2[3] = { "s_xx","s_yy","s_xy" };
  static const char* s3[6] = { "s_xx","s_yy","s_zz","s_xy","s_yz","s_xz" };

  static std::string label;
  if (prefix)
    label = prefix + std::string(" ");
  else
    label.clear();

  if (i >= (size_t)(nsd*(nsd+1)/2))
    label += "von Mises stress";
  else switch (nsd) {
  case 1: label += s1[i]; break;
  case 2: label += s2[i]; break;
  case 3: label += s3[i]; break;
  }

  return label.c_str();
}


NormBase* Elasticity::getNormIntegrand (TensorFunc* sol) const
{
  return new ElasticityNorm(*const_cast<Elasticity*>(this),sol);
}


bool ElasticityNorm::initElement (const std::vector<int>& MNPC)
{
  return problem.initElement(MNPC);
}


bool ElasticityNorm::initElementBou (const std::vector<int>& MNPC)
{
  return problem.initElementBou(MNPC);
}


size_t ElasticityNorm::getNoFields () const
{
  return anasol ? 4 : 2;
}


ElmNorm& ElasticityNorm::getElmNormBuffer (LocalIntegral*& elmInt)
{
  ElmNorm* eNorm = dynamic_cast<ElmNorm*>(elmInt);
  if (eNorm) return *eNorm;

  static double data[4];
  static ElmNorm buf(data);
  memset(data,0,4*sizeof(double));
  elmInt = &buf;
  return buf;
}


bool ElasticityNorm::evalInt (LocalIntegral*& elmInt, double detJW,
			      const Vector& N, const Matrix& dNdX,
			      const Vec3& X) const
{
  ElmNorm& pnorm = ElasticityNorm::getElmNormBuffer(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCmatrix(Cinv,X,true)) return false;

  // Evaluate the finite element stress field
  Vector sigmah, sigma;
  if (!problem.evalSol(sigmah,dNdX,X)) return false;

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += sigmah.dot(Cinv*sigmah)*detJW;

  if (problem.haveLoads())
  {
    // Evaluate the body load
    Vec3 f = problem.getBodyforce(X);
    // Evaluate the displacement field
    Vec3 u = problem.evalSol(N);
    // Integrate the external energy (f,u^h)
    pnorm[1] += f*u*detJW;
  }

  if (anasol)
  {
    // Evaluate the analytical stress field
    sigma = (*anasol)(X);
    // Integrate the energy norm a(u,u)
    pnorm[2] += sigma.dot(Cinv*sigma)*detJW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    sigma -= sigmah;
    pnorm[3] += sigma.dot(Cinv*sigma)*detJW;
  }

  return true;
}


bool ElasticityNorm::evalBou (LocalIntegral*& elmInt, double detJW,
			      const Vector& N, const Matrix&,
			      const Vec3& X, const Vec3& normal) const
{
  if (!problem.haveLoads()) return true;

  ElmNorm& pnorm = ElasticityNorm::getElmNormBuffer(elmInt);

  // Evaluate the surface traction
  Vec3 T = problem.getTraction(X,normal);
  // Evaluate the displacement field
  Vec3 u = problem.evalSol(N);

  // Integrate the external energy
  pnorm[1] += T*u*detJW;
  return true;
}
