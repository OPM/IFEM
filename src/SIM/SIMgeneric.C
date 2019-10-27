// $Id$
//==============================================================================
//!
//! \file SIMgeneric.C
//!
//! \date Aug 28 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generic SIM class with some added functionalities.
//!
//==============================================================================

#include "SIMgeneric.h"
#include "ModelGenerator.h"
#include "ASMbase.h"
#include "IntegrandBase.h"
#include "Utilities.h"
#include "IFEM.h"


ASMbase* SIMgeneric::createDefaultModel ()
{
  if (!myModel.empty()) return nullptr;

  ModelGenerator* gen = this->getModelGenerator(nullptr);
  if (!gen) return nullptr;

  bool okGen = gen->createGeometry(*this);
  nGlPatches = myModel.size();
  delete gen;

  return okGen && nGlPatches > 0 ? myModel.front() : nullptr;
}


Vector SIMgeneric::getSolution (const Vector& psol, const double* par,
                                int deriv, int patch) const
{
  if (psol.empty() || !par || opt.discretization < ASM::Spline)
    return Vector();

  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return Vector();

  size_t ndim = pch->getNoParamDim();
  std::vector<RealArray> params(ndim);
  for (size_t i = 0; i < ndim; i++)
    params[i].resize(1,par[i]);

  Matrix tmpVal;
  Vector localVec;
  pch->extractNodeVec(psol,localVec);
  if (!pch->evalSolution(tmpVal,localVec,&params.front(),false,deriv))
    return Vector();

  return tmpVal.getColumn(1);
}


int SIMgeneric::evalPoint (const double* xi, Vec3& X, double* param,
                           int patch, bool global) const
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return -1;

  double dummy[3] = {};
  int inod = pch->evalPoint(xi, param ? param : dummy, X);
  return inod > 0 && global ? pch->getNodeID(inod) : inod;
}


int SIMgeneric::findElementContaining (const double* param,
                                       int patch, bool global) const
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return -1;

  int iel = pch->findElementContaining(param);
  return iel > 0 && global ? pch->getElmID(iel) : iel;
}


double SIMgeneric::getReferenceNorm (const Vectors& gNorm, size_t adaptor) const
{
  if (gNorm.empty() || gNorm.front().empty())
    return 0.0;

  const Vector& fNorm = gNorm.front();
  if (adaptor < 1 && fNorm.size() > 2)
    return fNorm(3); // Using the analytical solution, |u|_ref = |u|
  else if (adaptor >= gNorm.size() || gNorm[adaptor].size() < 2)
    return -(double)adaptor; // Norm group index is out of range

  // |u|_ref = sqrt( |u^h|^2 + |e^*|^2 )
  return hypot(fNorm(1),gNorm[adaptor](2));
}


double SIMgeneric::getEffectivityIndex (const Vectors& gNorm,
                                        size_t idx, size_t inorm) const
{
  return gNorm[idx](inorm) / gNorm.front()(4);
}


size_t SIMgeneric::getVCPindex (size_t idx) const
{
  if (extrFunc.empty() || idx == 0)
    return 0;
  else if (idx > extrFunc.size())
    return 0;

  return (this->haveAnaSol() ? 4 : 2) + idx;
}


void SIMgeneric::printNorms (const Vectors& norms, size_t w) const
{
  if (norms.empty()) return;

  NormBase* norm = this->getNormIntegrand();
  const Vector& n = norms.front();

  IFEM::cout <<"Energy norm"
             << utl::adjustRight(w-11,norm->getName(1,1)) << n(1);
  if (n(2) != 0.0)
    IFEM::cout <<"\nExternal energy"
               << utl::adjustRight(w-15,norm->getName(1,2)) << n(2);

  if (this->haveAnaSol() && n.size() >= 4)
    IFEM::cout <<"\nExact norm"
               << utl::adjustRight(w-10,norm->getName(1,3)) << n(3)
               <<"\nExact error"
               << utl::adjustRight(w-11,norm->getName(1,4)) << n(4)
               <<"\nExact relative error (%) : "<< 100.0*n(4)/n(3);

  size_t i, j = 0, k = 0;
  while ((i = this->getVCPindex(++k)) && i <= n.size())
    IFEM::cout <<"\nVCP quantity"
               << utl::adjustRight(w-12,norm->getName(1,i)) << n(i);

  for (const SIMoptions::ProjectionMap::value_type& prj : opt.project)
    if (++j < norms.size())
      this->printNormGroup(norms[j],n,prj.second);

  IFEM::cout << std::endl;
  delete norm;
}
