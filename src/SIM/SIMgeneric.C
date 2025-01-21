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
#include "Function.h"
#include "SAM.h"
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


RealArray SIMgeneric::getSolution (const Vector& psol, const double* par,
                                   int deriv, int patch) const
{
  RealArray ptSol;
  if (psol.empty() || !par || opt.discretization < ASM::Spline)
    return ptSol;

  ASMbase* pch = this->getPatch(patch,true);
  if (!pch)
    return ptSol;

  size_t ndim = pch->getNoParamDim();
  std::vector<RealArray> params(ndim);
  for (size_t i = 0; i < ndim; i++)
    params[i].resize(1,par[i]);

  bool ok;
  Matrix tmpVal;
  Vector localVec;
  pch->extractNodeVec(psol,localVec);
  if (myProblem->getIntegrandType() & Integrand::PIOLA_MAPPING)
    ok = pch->evalSolutionPiola(tmpVal,localVec,&params.front(),false);
  else
    ok = pch->evalSolution(tmpVal,localVec,&params.front(),false,deriv);
  if (ok)
    ptSol = tmpVal.getColumn(1);

  return ptSol;
}


RealArray SIMgeneric::getInterfaceForces (const RealArray& sf,
                                          const RealArray& weights,
                                          int code) const
{
  RealArray force(nsd,0.0);
  if (!mySam)
    return force;

  IntVec glbNodes;
  this->getBoundaryNodes(code,glbNodes);
#if SP_DEBUG > 1
  std::cout <<"\nInternal nodal forces at interface "<< code
            <<" with "<< glbNodes.size() <<" nodes"<< std::endl;
#endif

  for (int inod : glbNodes)
  {
    double w = inod <= (int)weights.size() ? weights[inod-1] : 1.0;
    std::pair<int,int> dof = mySam->getNodeDOFs(inod);
#if SP_DEBUG > 1
    std::cout <<"Node "<< inod <<":";
    for (int j = dof.first-1; j < dof.second; j++)
      std::cout <<" "<< sf[j];
    std::cout << std::endl;
#endif
    for (unsigned char i = 0; i < nsd; i++, dof.first++)
      if (dof.first <= dof.second && dof.first < (int)sf.size())
        force[i] += w*sf[dof.first-1];
      else
        break;
  }

  return force;
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
  else if (!this->haveAnaSol() && idx > extrFunc.size())
    return 0;
  else if (idx > 2*extrFunc.size())
    return 0;

  size_t iVol = this->getVolumeIndex();
  // Assume here that the VCP quantities follow directly after the volume/area
  return iVol ? iVol + idx : 0;
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

  if (norms.size() > 1+opt.project.size())
  {
    if (norms.back().size() == 1)
      IFEM::cout <<"\nVol(D) : "<<  norms.back()[0];
    else for (i = 0; i < norms.back().size(); i++)
      IFEM::cout <<"\nVol(D"<< i+1 <<") : "<<  norms.back()[i];
  }

  for (const SIMoptions::ProjectionMap::value_type& prj : opt.project)
    if (++j < norms.size())
      this->printNormGroup(norms[j],n,prj.second);

  IFEM::cout << std::endl;
  delete norm;
}


/*!
  This method is used to fix the problem that all norm quantities (global and
  on element level) will be subjected to a final square-root operation after
  they have been integrated and accumulated, through NormBase::applyFinalOp().
  However, a few quantities (like the volume and variationally consistent
  postprocessing values) should not reveive this operation.
*/

bool SIMgeneric::revertSqrt (Vectors& gNorm, Matrix* eNorm)
{
  if (gNorm.empty())
    return false;

  size_t iVol = this->getVolumeIndex();
  if (iVol > gNorm.front().size())
    iVol = 0;
  else if (iVol > 0)
  {
    // Undo the square-root final operation for the global and element volumes
    gNorm.front()(iVol) *= gNorm.front()(iVol);
    if (!eNorm || iVol > eNorm->rows())
      iVol = 0;
    else for (size_t j = 1; j <= eNorm->cols(); j++)
      (*eNorm)(iVol,j) *= (*eNorm)(iVol,j);
  }

  // Check for variationally consistent postprocessing (VCP) quantities
  size_t i, k = 0;
  while ((i = this->getVCPindex(++k)))
  {
    double Dvol = 0.0; // Volume of support of the extraction function
    if (eNorm && i <= eNorm->rows())
      for (size_t j = 1; j <= eNorm->cols(); j++)
      {
        double& vcpq = (*eNorm)(i,j);
        // Check if element is within the support of the extraction function
        if (fabs(vcpq) > 1.0e-12 && iVol > 0)
          Dvol += (*eNorm)(iVol,j);

        // Undo the square-root final operation for the element quantity,
        vcpq = copysign(vcpq*vcpq,vcpq); // while preserving the sign
      }

    // Store the domain volume in gNorm
    size_t j = this->haveAnaSol() ? (k-1)/2 : k-1;
    if (j < gNorm.back().size())
      gNorm.back()[j] = Dvol*Dvol;

    if (i <= gNorm.front().size())
    {
      double& vcpq = gNorm.front()(i);
      // Divide by the volume to obtain volume-averaged quantity
      // if we are doing point-wise result extraction.
      // Notice that this is done _before_ the square-ing.
      // This is because the final sqrt-operation has not been
      // performed yet on the global norm quantities.
      if (Dvol > 0.0 && extrFunc[j]->getType() == 3)
        vcpq /= Dvol;

      // Undo the square-root final operation for the global quantity,
      vcpq = copysign(vcpq*vcpq,vcpq); // while preserving the sign
    }
  }

  return true;
}
