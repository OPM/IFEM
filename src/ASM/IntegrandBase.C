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
#include "ASMbase.h"
#include "GlobalIntegral.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TensorFunction.h"
#include "Utilities.h"
#include "Field.h"
#include "Fields.h"
#include <cstring>
#include <cstdio>


/*!
  The default implementation of this method allocates for one solution vector
  in \ref primsol when the \a mode is set to SIM::RECOVERY or SIM::NORMS,
  but only if the vector is empty.
  This is sufficient for most linear problems with a single right-hand side.
*/

void IntegrandBase::setMode (SIM::SolutionMode mode)
{
  if ((m_mode = mode) >= SIM::RECOVERY && primsol.empty())
    primsol.resize(1);
}


/*!
  Override this method if the integrand needs some patch-specific data
  to be initialized before performing the numerical integration.
  The default version only passes the patch index to the initPatch() method.
*/

void IntegrandBase::initForPatch (const ASMbase* pch)
{
  if (pch) this->initPatch(pch->idx);
}


/*!
  Override this method if the problem to be solved consists of multiple
  IntegrandBase objects. This method is then supposed to return the
  corresponding integrated quantity for this integrand.
  The default implementation returns the object provided as argument.
*/

GlobalIntegral& IntegrandBase::getGlobalInt (GlobalIntegral* gq) const
{
  if (gq) return *gq;

  std::cerr <<"  ** IntegrandBase::getGlobalInt: No global integral"
            <<" - using dummy."<< std::endl;

  static GlobalIntegral dummy;
  return dummy;
}


/*!
  The default implementation returns an ElmMats object with one left-hand-side
  matrix (unless we are doing a boundary integral) and one right-hand-side
  vector. The dimension of the element matrices are assumed to be \a npv*nen.
  Override this method if your integrand needs more element matrices.
*/

LocalIntegral* IntegrandBase::getLocalIntegral (size_t nen, size_t,
                                                bool neumann) const
{
  ElmMats* result = new ElmMats(!neumann && m_mode < SIM::RECOVERY);
  result->rhsOnly = m_mode >= SIM::RHS_ONLY;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(npv*nen);

  return result;
}


/*!
  This method can be used when only the first primary solution vector is needed.
  Thus, the returned \a elmVec array of element solution vectors will always
  have dimension 1.

  The nodal point correspondance array \a MNPC may contain more nodal indices
  than should be used in the extraction process. In that case, the number of
  nodes to omit at the end needs to be specified through the \a nskip argument.
*/

bool IntegrandBase::initElement1 (const std::vector<int>& MNPC,
                                  Vectors& elmVec, size_t nskip) const
{
  elmVec.resize(1);
  if (primsol.empty() || primsol.front().empty())
    return true; // No solution fields yet, return an empty vector

  // Extract the first primary solution vector for this element
  int ierr = utl::gather(MNPC,npv,primsol.front(),elmVec.front(),0,nskip);
  if (ierr > 0)
  {
    std::cerr <<" *** IntegrandBase::initElement1: Detected "
              << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

#if SP_DEBUG > 2
  std::cout <<"Element solution vector"<< elmVec.front();
#endif
  return true;
}


/*!
  This is the basic implementation of the element solution vector extraction.
  The output array \a elmVec will contain one element-level solution vector
  corresponding to each patch-level vector in \ref primsol.

  The nodal point correspondance array \a MNPC may contain more nodal indices
  than should be used in the extraction process. In that case, the number of
  nodes to omit at the end needs to be specified through the \a nskip argument.
*/

bool IntegrandBase::initElement2 (const std::vector<int>& MNPC,
                                  Vectors& elmVec, size_t nskip) const
{
  size_t nsol = primsol.size();
  while (nsol > 1 && primsol[nsol-1].empty()) nsol--;
  if (nsol <= 1) // Only one (or none) primary solution vector to extract
    return this->initElement1(MNPC,elmVec,nskip);

  // Extract all primary solution vectors for this element
  int ierr = 0;
  elmVec.resize(nsol);
  for (size_t i = 0; i < nsol && ierr == 0; i++)
    if (!primsol[i].empty())
      ierr = utl::gather(MNPC,npv,primsol[i],elmVec[i],0,nskip);

#if SP_DEBUG > 2
  for (size_t j = 0; j < nsol; j++)
    std::cout <<"Element solution vector "<< j+1 << elmVec[j];
#endif

  if (ierr == 0) return true;

  std::cerr <<" *** IntegrandBase::initElement2: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


/*!
  The default implementation extracts the element-level solution vectors,
  stored in the provided \a elmInt argument, corresponding to each of the
  patch-wise solution vectors in the member \ref primsol.
  Override this method if your integrand requires some additional or other
  element initializations.
*/

bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 LocalIntegral& elmInt)
{
  return this->initElement2(MNPC,elmInt.vec);
}


/*!
  Override this method if your integrand needs either the FiniteElement object
  itself, the element center, or the total number of integration points within
  the element. The default implementation forwards to the overloaded method
  not taking any of \a fe, \a X0 or \a nPt as arguments.
*/

bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 const FiniteElement&, const Vec3&, size_t,
                                 LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
}


/*!
  This method is used when a mixed interpolation basis is used. The default
  implementation forwards to the single-basis version, using first basis only.
  Override this method for mixed problems requiring element-level solution
  vectors also for the other bases.
*/

bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 const std::vector<size_t>& elem_sizes,
                                 const std::vector<size_t>& basis_sizes,
                                 LocalIntegral& elmInt)
{
  std::vector<int> MNPC1(MNPC.begin(), MNPC.begin()+elem_sizes.front());
  return this->initElement(MNPC1,elmInt);
}


bool IntegrandBase::initElement (const std::vector<int>& MNPC,
                                 const MxFiniteElement&,
                                 const std::vector<size_t>& elem_sizes,
                                 const std::vector<size_t>& basis_sizes,
                                 LocalIntegral& elmInt)
{
  return this->initElement(MNPC, elem_sizes, basis_sizes, elmInt);
}


/*!
  The default implementation extracts the element-level vector only for the
  first (current) primary solution vector.
  Override this method if your boundary integrand requires more than that.
*/

bool IntegrandBase::initElementBou (const std::vector<int>& MNPC,
                                    LocalIntegral& elmInt)
{
  return this->initElement1(MNPC,elmInt.vec);
}


/*!
  This method is used when a mixed interpolation basis is used. The default
  implementation forwards to the single-basis version, using first basis only.
  Override this method for mixed problems requiring element-level solution
  vectors also for the other bases.
*/

bool IntegrandBase::initElementBou (const std::vector<int>& MNPC,
                                    const std::vector<size_t>& elem_sizes,
                                    const std::vector<size_t>& basis_sizes,
                                    LocalIntegral& elmInt)
{
  std::vector<int> MNPC1(MNPC.begin(),MNPC.begin()+elem_sizes.front());
  return this->initElementBou(MNPC1,elmInt);
}


/*!
  This method can be used when only the first primary solution vector is needed
  in the secondary solution evaluation. Thus, only one element-level solution
  vector is extracted and forwarded to the virtual method evalSol2().
*/

bool IntegrandBase::evalSol1 (Vector& s, const FiniteElement& fe, const Vec3& X,
                              const std::vector<int>& MNPC, size_t nskip) const
{
  if (primsol.empty() || primsol.front().empty())
  {
    std::cerr <<" *** IntegrandBase::evalSol: No solution vector."<< std::endl;
    return false;
  }

  // Extract the first primary solution vector for this element
  Vectors elmVec(1);
  int ierr = utl::gather(MNPC,npv,primsol.front(),elmVec.front(),0,nskip);
  if (ierr > 0)
  {
    std::cerr <<" *** IntegrandBase::evalSol: Detected "
              << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

  return this->evalSol2(s,elmVec,fe,X);
}


/*!
  This method must be overridded in a sub-class.
*/

bool IntegrandBase::evalSol2 (Vector&, const Vectors&,
                              const FiniteElement&, const Vec3&) const
{
  std::cerr <<" *** IntegrandBase::evalSol2: Not implemented."<< std::endl;
  return false;
}


/*!
  The default implementation assumes that only the first primary solution
  vector is needed in the secondary solution evaluation, and thus forwards
  to the evalSol1() method.
  Override this method if the secondary solution depends on more than that.
*/

bool IntegrandBase::evalSol (Vector& s, const FiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC) const
{
  return this->evalSol1(s,fe,X,MNPC);
}


/*!
  This method is used when a mixed interpolation basis is used. The default
  implementation forwards to the single-basis version, using first basis only.
  Override this method for mixed problems requiring element-level solution
  vectors also for the other bases in the secondary solution evaluation.
*/

bool IntegrandBase::evalSol (Vector& s, const MxFiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC,
                             const std::vector<size_t>& elem_sizes,
                             const std::vector<size_t>& basis_sizes) const
{
  std::vector<int> MNPC1(MNPC.begin(),MNPC.begin()+elem_sizes.front());
  return this->evalSol(s,fe,X,MNPC1);
}


bool IntegrandBase::evalSol (Vector& s, const TensorFunc& asol,
                             const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd*nsd);
  return true;
}


bool IntegrandBase::evalSol (Vector& s, const STensorFunc& asol,
                             const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd*(nsd+1)/2);
  return true;
}


bool IntegrandBase::evalSol (Vector& s, const VecFunc& asol,
                             const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


void IntegrandBase::resetSolution ()
{
  for (Vector& sol : primsol) sol.clear();
}


void IntegrandBase::printSolution (std::ostream& os, int pindx)
{
  if (std::all_of(primsol.begin(), primsol.end(),
                  [](const Vector& sol) { return sol.empty(); }))
    return; // No solution, or all-empty

  int isol = 0;
  os <<"\nCurrent solution for Patch "<< pindx;
  for (Vector& sol : primsol)
  {
    isol++;
    if (sol.empty()) continue;

    if (primsol.size() > 1)
      os <<"\nSolution vector "<< isol;

    if (npv > 0)
    {
      int idof = 0, node = 0;
      for (double sval : sol)
        if (npv < 2 || (++idof)%npv == 1)
          os <<"\nNode "<< ++node <<": "<< sval;
        else
          os <<" "<< sval;
    }
    else // For mixed problems there is a variable number of DOFs per node
      os << sol; // so just print the raw solution vector instead

    os << std::endl;
  }
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


std::string IntegrandBase::getField1Name (size_t idx, const char* prefix) const
{
  char name[32];
  sprintf(name,"primary solution %zu",1+idx);
  return prefix ? prefix +  std::string(" ") + name : name;
}


std::string IntegrandBase::getField2Name (size_t idx, const char* prefix) const
{
  char name[32];
  sprintf(name,"secondary solution %zu",1+idx);
  return prefix ? prefix + std::string(" ") + name : name;
}


bool IntegrandBase::inActive (int iel) const
{
  if (elmGrp.empty()) return false;

  return std::find(elmGrp.begin(),elmGrp.end(),iel) == elmGrp.end();
}


NormBase::~NormBase ()
{
  for (Fields* f : prjFld)
    delete f;
}


void NormBase::initProjection (size_t nproj)
{
  prjFld.resize(nproj,nullptr);
  prjsol.resize(nproj);
}


Vector& NormBase::getProjection (size_t i)
{
  if (i < prjsol.size())
    return prjsol[i];

  static Vector dummy;
  return dummy;
}


void NormBase::setProjectedFields (Fields* f, size_t idx)
{
  if (idx < prjFld.size())
  {
    std::swap(f,prjFld[idx]);
    delete f;
  }
}


bool NormBase::initProjection (const std::vector<int>& MNPC,
                               LocalIntegral& elmInt, size_t nExtraNodes)
{
  // Extract projected solution vectors for this element
  Vectors& psol = static_cast<ElmNorm&>(elmInt).psol;
  psol.resize(prjsol.size());

  int ierr = 0;
  for (size_t i = 0; i < prjsol.size() && ierr == 0; i++)
    if (!prjsol[i].empty())
      ierr = utl::gather(MNPC,nrcmp,prjsol[i],psol[i],0,nExtraNodes);

  if (ierr == 0) return true;

  std::cerr <<" *** NormBase::initProjection: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;
  return false;
}


LocalIntegral* NormBase::getLocalIntegral (size_t, size_t iEl, bool) const
{
  if (myProblem.inActive(iEl)) return nullptr;

  if (lints && iEl > 0 && iEl <= lints->size()) return (*lints)[iEl-1];

  // Element norms are not requested, so allocate one internally instead,
  // that will delete itself when invoking the destruct method.
  size_t j, norms = 0;
  size_t groups = this->getNoFields(0);
  for (j = 1; j <= groups; j++)
    norms += this->getNoFields(j);

  return new ElmNorm(norms);
}


double NormBase::applyFinalOp (double value) const
{
  switch (finalOp)
    {
    case ASM::ABS:
      return fabs(value);
    case ASM::SQRT:
      return value < 0.0 ? -sqrt(-value) : sqrt(value);
    default:
      return value;
    }
}


void NormBase::addBoundaryTerms (Vectors& gNorm, double energy) const
{
  if (gNorm.empty() || gNorm.front().size() < 2 || energy == 0.0) return;

  double& extEnergy = gNorm.front()[1];
#ifdef SP_DEBUG
  if (extEnergy != 0.0)
    std::cout <<"External energy contribution from interior terms: "
              << this->applyFinalOp(extEnergy) << std::endl;
  std::cout <<"External energy contribution from boundary terms: "
            << this->applyFinalOp(energy) << std::endl;
#endif
  extEnergy += energy;
}


size_t NormBase::getNoFields (int group) const
{
  return group > 0 ? 0 : 1 + prjsol.size();
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            const FiniteElement& fe,
                            const Vec3& Xc, size_t nPt,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt,myProblem.getNoGLMs()) &&
         myProblem.initElement(MNPC,fe,Xc,nPt,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt,myProblem.getNoGLMs()) &&
         myProblem.initElement(MNPC,elmInt);
}


bool NormBase::initElement (const std::vector<int>& MNPC,
                            const std::vector<size_t>& elem_sizes,
                            const std::vector<size_t>& basis_sizes,
                            LocalIntegral& elmInt)
{
  size_t nExtraNodes = MNPC.size() - elem_sizes.front();
  return this->initProjection(MNPC,elmInt,nExtraNodes) &&
         myProblem.initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC,
                               LocalIntegral& elmInt)
{
  if (projBou && !this->initProjection(MNPC,elmInt,myProblem.getNoGLMs()))
    return false;

  return myProblem.initElementBou(MNPC,elmInt);
}


bool NormBase::initElementBou (const std::vector<int>& MNPC,
                               const std::vector<size_t>& elem_sizes,
                               const std::vector<size_t>& basis_sizes,
                               LocalIntegral& elmInt)
{
  size_t nExtraNodes = MNPC.size() - elem_sizes.front();
  if (projBou && !this->initProjection(MNPC,elmInt,nExtraNodes))
    return false;

  return myProblem.initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


std::string NormBase::getName (size_t i, size_t j, const char* prefix) const
{
  char comp[32];
  sprintf(comp,"norm_%zu.%zu",i,j);
  if (!prefix) return comp;

  return prefix + std::string(" ") + comp;
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
  this->clearBuffer();
}


void ForceBase::clearBuffer ()
{
  for (LocalIntegral* lint : eForce)
    delete lint;
  eForce.clear();
  delete[] eBuffer;
  eBuffer = nullptr;
}


void ForceBase::initBuffer (size_t nel)
{
  size_t ncmp = this->getNoComps();

  if (eForce.size() != nel) {
    this->clearBuffer();

    eBuffer = new double[ncmp*nel];

    eForce.reserve(nel);
    double* q = eBuffer;
    for (size_t i = 0; i < nel; i++, q += ncmp)
      eForce.push_back(new ElmNorm(q,ncmp));
  }

  memset(eBuffer,0,ncmp*nel*sizeof(double));
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
  for (LocalIntegral* lint : eForce)
  {
    ElmNorm& elmf = *static_cast<ElmNorm*>(lint);
    for (size_t j = 0; j < gForce.size() && j < elmf.size(); j++)
      gForce[j] += elmf[j];
  }
}
