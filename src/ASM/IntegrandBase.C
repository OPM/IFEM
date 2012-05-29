// $Id$
//==============================================================================
//!
//! \file IntegrandBase.C
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base classes representing FEM integrands.
//!
//==============================================================================

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "GlbNorm.h"
#include "Utilities.h"


void IntegrandBase::resetSolution ()
{
  for (size_t i = 0; i < primsol.size(); i++)
    primsol[i].clear();
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC,
				 const Vec3&, size_t, LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 LocalIntegral& elmInt)
{
  // Extract the primary solution vectors for this element
  int ierr = 0;
  elmInt.vec.resize(primsol.empty() ? 1 : primsol.size());
  for (size_t i = 0; i < primsol.size() && ierr == 0; i++)
    if (!primsol[i].empty())
      ierr = utl::gather(MNPC,npv,primsol[i],elmInt.vec[i]);

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool IntegrandBase::initElement (const std::vector<int>&,
				 const std::vector<int>&, size_t,
                                 LocalIntegral&)
{
  std::cerr <<" *** IntegrandBase::initElement must be reimplemented"
	    <<" for mixed problems."<< std::endl;
  return false;
}


bool IntegrandBase::initElementBou (const std::vector<int>& MNPC,
                                    LocalIntegral& elmInt)
{
  // Extract current primary solution vector for this element
  int ierr = 0;
  elmInt.vec.resize(1);
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),elmInt.vec.front());

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElementBou: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool IntegrandBase::initElementBou (const std::vector<int>&,
				    const std::vector<int>&, size_t,
                                    LocalIntegral&)
{
  std::cerr <<" *** IntegrandBase::initElementBou must be reimplemented"
	    <<" for mixed problems."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const Vector&, const Matrix&,
			     const Vec3&, const std::vector<int>&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector& s, const Vector& N,
			     const Matrix& dNdX, const Matrix3D&,
			     const Vec3& X, const std::vector<int>& MNPC) const
{
  return this->evalSol(s,N,dNdX,X,MNPC);
}


bool IntegrandBase::evalSol (Vector& s, const Vector& N1, const Vector&,
			     const Matrix& dN1dX, const Matrix&, const Vec3& X,
			     const std::vector<int>& MNPC1,
			     const std::vector<int>&) const
{
  return this->evalSol(s,N1,dN1dX,X,MNPC1);
}


bool IntegrandBase::evalSol (Vector&, const TensorFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const STensorFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


bool IntegrandBase::evalSol (Vector&, const VecFunc&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol not implemented."<< std::endl;
  return false;
}


void IntegrandBase::registerVector (const std::string& name, Vector* vec)
{
  myFields[name] = vec;
}


Vector* IntegrandBase::getNamedVector (const std::string& name) const
{
  std::map<std::string,Vector*>::const_iterator it = myFields.find(name);
  return it == myFields.end() ? NULL : it->second;
}


void NormBase::initIntegration (const TimeDomain& time)
{
  myProblem.initIntegration(time);
}


Vector& NormBase::getProjection (size_t i)
{
  if (i < prjsol.size())
    return prjsol[i];

  static Vector dummy;
  return dummy;
}


bool NormBase::initProjection (const std::vector<int>& MNPC,
                               LocalIntegral& elmInt)
{
  // Extract projected solution vectors for this element
  Vectors& psol = static_cast<ElmNorm&>(elmInt).psol;
  psol.resize(prjsol.size());

  int ierr = 0;
  for (size_t i = 0; i < prjsol.size() && ierr == 0; i++)
    if (!prjsol[i].empty())
      ierr = utl::gather(MNPC,nrcmp,prjsol[i],psol[i]);

  if (ierr == 0) return true;

  std::cerr <<" *** NormBase::initProjection: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;
  return false;
}


LocalIntegral* NormBase::getLocalIntegral (size_t, size_t iEl, bool) const
{
  if (lints && iEl > 0 && iEl <= lints->size()) return (*lints)[iEl-1];

  // Element norms are not requested, so allocate one internally instead that
  // will delete itself when invoking the destruct method.
  return new ElmNorm(this->getNoFields());
}


bool NormBase::initElement (const std::vector<int>& MNPC,
			    const Vec3& Xc, size_t nPt,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,Xc,nPt,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC1,
			    const std::vector<int>& MNPC2, size_t n1,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC1,elmInt) &&
         myProblem.initElement(MNPC1,MNPC2,n1,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC,
                               LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC1,
			       const std::vector<int>& MNPC2, size_t n1,
                               LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC1,MNPC2,n1,elmInt);
}


const char* NormBase::getName (size_t& i, bool haveExact, const char* prefix)
{
  static const char* s[7] = {
    "a(u^h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "a(u^r,u^r)^0.5",
    "a(e',e')^0.5, e'=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "(a(e',e')/a(e,e))^0.5"
  };

  if (!haveExact) // Skip norms present only together with exact solutions
    if (i == 1)
      i = 3;
    else if (i == 5)
      i = 7;

  if (!prefix && i < 7) return s[i];

  static std::string name;
  if (prefix)
    name = prefix + std::string(" ");
  else
    name.clear();

  if (i < 7)
    name += s[i];
  else
  {
    char comp[16];
    sprintf(comp,"norm_%lu",1+i);
    name += comp;
  }

  return name.c_str();
}


ForceBase::~ForceBase ()
{
  for (size_t i = 0; i < eForce.size(); i++)
    delete eForce[i];
  delete[] eBuffer;
}


bool ForceBase::initBuffer (size_t nel)
{
  if (eBuffer) return false;

  size_t ncmp = this->getNoComps();
  eBuffer = new double[ncmp*nel];
  memset(eBuffer,0,ncmp*nel*sizeof(double));

  eForce.reserve(nel);
  double* q = eBuffer;
  for (size_t i = 0; i < nel; i++, q += ncmp)
    eForce.push_back(new ElmNorm(q,ncmp));

  return true;
}


LocalIntegral* ForceBase::getLocalIntegral (size_t, size_t iEl, bool) const
{
  if (iEl > 0 && iEl <= eForce.size()) return eForce[iEl-1];

  // No internal buffers. Allocate an ElmNorm object that will delete itself
  // when invoking its destruct method.
  return new ElmNorm(this->getNoComps());
}


void ForceBase::initIntegration (const TimeDomain& time)
{
  myProblem.initIntegration(time);
}


bool ForceBase::initElementBou (const std::vector<int>& MNPC,
                                LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,elmInt);
}


bool ForceBase::initElementBou (const std::vector<int>& MNPC1,
                                const std::vector<int>& MNPC2, size_t n1,
                                LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC1,MNPC2,n1,elmInt);
}


void ForceBase::assemble (GlbNorm& force) const
{
  for (size_t i = 0; i < eForce.size(); i++)
    force.assemble(eForce[i]);
}
