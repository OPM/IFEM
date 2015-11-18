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
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Utilities.h"
#include "Field.h"
#include "Fields.h"
#include <cstring>
#include <sstream>
#include <cstdio>


/*!
  \brief Convenience function writing error message for non-implemented methods.
*/

static bool Ierror (const char* name)
{
  std::cerr <<" *** IntegrandBase::"<< name <<" not implemented."<< std::endl;
  return false;
}


/*!
  The default implementation returns an ElmMats object with one left-hand-side
  matrix (unless we are doing a boundary integral) and one right-hand-side
  vector. The dimension of the element matrices are assumed to be \a npv*nen.
  Reimplement this method if your integrand needs more element matrices.
*/

LocalIntegral* IntegrandBase::getLocalIntegral (size_t nen, size_t,
                                                bool neumann) const
{
  ElmMats* result = new ElmMats(!neumann && m_mode != SIM::RECOVERY);
  result->rhsOnly = m_mode == SIM::RHS_ONLY || m_mode == SIM::RECOVERY;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(npv*nen);

  return result;
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 LocalIntegral& elmInt)
{
  // Extract all primary solution vectors for this element
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


/*!
  Reimplement this method if your integrand need either the FiniteElement object
  itself, the element center, or the total number of integration points within
  the element. The default implementation forwards to an overloaded method
  not taking any of \a fe, \a X0 or \a nPt as arguments.
*/

bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 const FiniteElement&, const Vec3&, size_t,
                                 LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
}


/*!
  The default implementation forwards to the single-basis version.
*/

bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 const std::vector<size_t>& elem_sizes,
                                 const std::vector<size_t>& basis_sizes,
                                 LocalIntegral& elmInt)
{
  std::vector<int> MNPC1(MNPC.begin(), MNPC.begin()+elem_sizes.front());
  return this->initElement(MNPC1,elmInt);
}


bool IntegrandBase::initElementBou (const std::vector<int>& MNPC,
                                    LocalIntegral& elmInt)
{
  // Extract (only) the current primary solution vector for this element
  int ierr = 0;
  elmInt.vec.resize(1);
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),elmInt.vec.front());

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElementBou: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


/*!
  The default implementation forwards to the single-basis version.
*/

bool IntegrandBase::initElementBou (const std::vector<int>& MNPC,
                                    const std::vector<size_t>& elem_sizes,
                                    const std::vector<size_t>& basis_sizes,
                                    LocalIntegral& elmInt)
{
  std::vector<int> MNPC1(MNPC.begin(), MNPC.begin()+elem_sizes.front());
  return this->initElementBou(MNPC1,elmInt);
}


bool IntegrandBase::evalSol (Vector& s, const FiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC) const
{
  return Ierror("evalSol(Vector&,const FiniteElement& fe,const Vec3&,...)");
}


bool IntegrandBase::evalSol (Vector& s, const MxFiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC,
                             const std::vector<size_t>& elem_sizes) const
{
  std::vector<int> MNPC1(MNPC.begin(), MNPC.begin()+elem_sizes.front());
  return this->evalSol(s,fe,X,MNPC1);
}


bool IntegrandBase::evalSol (Vector& s, const TensorFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


bool IntegrandBase::evalSol (Vector& s, const STensorFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


bool IntegrandBase::evalSol (Vector& s, const VecFunc& asol, const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


void IntegrandBase::resetSolution ()
{
  for (Vectors::iterator it = primsol.begin(); it != primsol.end(); ++it)
    it->clear();
}


void IntegrandBase::registerVector (const std::string& name, Vector* vec)
{
  myFields[name] = vec;
}


void IntegrandBase::setNamedField (const std::string&, Field* f)
{
  delete f;
}


void IntegrandBase::setNamedFields (const std::string&, Fields* f)
{
  delete f;
}


Vector* IntegrandBase::getNamedVector (const std::string& name) const
{
  std::map<std::string,Vector*>::const_iterator it = myFields.find(name);
  return it == myFields.end() ? nullptr : it->second;
}

std::string IntegrandBase::getField1Name(size_t idx, const char* prefix) const 
{
  std::stringstream stream;
  if(prefix)
    stream << prefix << " ";
  stream << "primary solution " << idx+1;
  return stream.str();
}

std::string IntegrandBase::getField2Name(size_t idx, const char* prefix) const 
{
  std::stringstream stream;
  if(prefix)
    stream << prefix << " ";
  stream << "secondary solution " << idx+1;
  return stream.str();
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
  size_t norms = 0;
  size_t groups = this->getNoFields(0);
  for (size_t j = 1; j <= groups; ++j)
    norms += this->getNoFields(j);

  return new ElmNorm(norms);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            const FiniteElement& fe,
                            const Vec3& Xc, size_t nPt,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,fe,Xc,nPt,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            const std::vector<size_t>& elem_sizes,
                            const std::vector<size_t>& basis_sizes,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC,
                               LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC,
                               const std::vector<size_t>& elem_sizes,
                               const std::vector<size_t>& basis_sizes,
                               LocalIntegral& elmInt)
{
  return myProblem.initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


std::string NormBase::getName (size_t i, size_t j, const char* prefix) const
{
  char comp[32];
  sprintf(comp,"norm_%lu.%lu",i,j);
  if (prefix)
  {
    std::string tmp(comp);
    if (strlen(prefix)+1+tmp.size() < sizeof(comp))
      strcpy(comp,prefix);
    else
      strncpy(comp,prefix,sizeof(comp)-1-tmp.size());
    strcat(comp," ");
    strcat(comp,tmp.c_str());
  }

  return comp;
}


int NormBase::getIntegrandType () const
{
  // Mask off the element interface flag, if set
  return myProblem.getIntegrandType() & ~INTERFACE_TERMS;
}


int NormBase::getReducedIntegration (int n) const
{
  return myProblem.getReducedIntegration(n);
}


bool NormBase::reducedInt (LocalIntegral& elmInt,
                           const FiniteElement& fe, const Vec3& X) const
{
  return myProblem.reducedInt(elmInt,fe,X);
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


bool ForceBase::initElementBou (const std::vector<int>& MNPC,
                                LocalIntegral& elmInt)
{
  // Note that we invoke initElement (and not initElementBou) of the problem
  // integrand here, because the forces may depend on all solution variables.
  return myProblem.initElement(MNPC,elmInt);
}


bool ForceBase::initElementBou (const std::vector<int>& MNPC,
                                const std::vector<size_t>& elem_sizes,
                                const std::vector<size_t>& basis_sizes,
                                LocalIntegral& elmInt)
{
  // Note that we invoke initElement (and not initElementBou) of the problem
  // integrand here, because the forces may depend on all solution variables.
  return myProblem.initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


void ForceBase::assemble (RealArray& gForce) const
{
  gForce.resize(this->getNoComps());
  std::fill(gForce.begin(),gForce.end(),0.0);
  for (size_t i = 0; i < eForce.size(); i++)
  {
    ElmNorm& elmf = *static_cast<ElmNorm*>(eForce[i]);
    for (size_t j = 0; j < gForce.size() && j < elmf.size(); j++)
      gForce[j] += elmf[j];
  }
}
