// $Id$
//==============================================================================
//!
//! \file ForceIntegrator.C
//!
//! \date May 10 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for integration of boundary and nodal forces.
//!
//==============================================================================

#include "ForceIntegrator.h"
#include "SIMbase.h"
#include "ASMbase.h"
#include "IntegrandBase.h"
#include "GlbForceVec.h"
#include "Vec3.h"


Vector SIM::getBoundaryForce (const Vectors& solution, SIMbase* model, int code,
                              const TimeDomain& time, const Vector* X0)
{
  if (X0)
  {
    Vec3 X(X0->ptr(),X0->size());
    return SIM::getBoundaryForce(solution,model,code,time,&X);
  }
  return SIM::getBoundaryForce(solution,model,code,time);
}


Vector SIM::getBoundaryForce (const Vectors& solution, SIMbase* model, int code,
                              const TimeDomain& time, const Vec3* X0)
{
  ForceBase* forceInt = model->getBoundaryForceIntegrand(X0);
  if (!forceInt)
  {
    std::cerr <<" *** SIM::getBoundaryForce: No force integrand."<< std::endl;
    return Vector();
  }

  forceInt->initBuffer(model->getNoElms());

  // Integrate forces for given boundary segment
  bool ok = true;
  size_t prevPatch = 0;
  GlobalIntegral dummy;
  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop() && ok; p++)
    if (abs(p->pindx) == code)
    {
      ASMbase* patch = model->getPatch(p->patch);
      if (!patch)
        ok = false;
      else if (abs(p->ldim)+1 == patch->getNoParamDim())
      {
        if (p->patch != prevPatch)
          ok = model->extractPatchSolution(solution,p->patch-1);
        ok &= patch->integrate(*forceInt,abs(p->lindx),dummy,time);
        prevPatch = p->patch;
      }
    }

  // Assemble the element force contributions into the global force resultant
  Vector force;
  forceInt->assemble(force);
  delete forceInt;

  if (ok) return force;

  std::cerr <<" *** SIM::getBoundaryForce: Failed to evaluate boundary force"
            << std::endl;
  return Vector();
}


bool SIM::getNodalForces (const Vectors& solution, SIMbase* model, int code,
                          const TimeDomain& time, GlbForceVec& force)
{
  force.initialize();

  ForceBase* forceInt = model->getNodalForceIntegrand();
  if (!forceInt)
  {
    std::cerr <<" *** SIM::getNodalForces: No force integrand."<< std::endl;
    return false;
  }

  // Integrate nodal forces for given boundary segment
  bool ok = true;
  size_t prevPatch = 0;
  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop() && ok; p++)
    if (abs(p->pindx) == code)
    {
      ASMbase* patch = model->getPatch(p->patch);
      if (!patch)
        ok = false;
      else if (abs(p->ldim)+1 == patch->getNoParamDim())
      {
        if (p->patch != prevPatch)
          ok = model->extractPatchSolution(solution,p->patch-1);
        ok &= patch->integrate(*forceInt,abs(p->lindx),force,time);
        prevPatch = p->patch;
      }
    }

  delete forceInt;
  if (ok & force.finalize())
    return true;

  std::cerr <<" *** SIM::getNodalForces: Failed to evaluate nodal forces"
            << std::endl;
  return false;
}


bool SIM::initBoundaryNodeMap (SIMbase* model, int code, GlbForceVec& force)
{
  ASMbase* patch;
  IntVec glbNodes;
  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop(); p++)
    if (abs(p->pindx) == code && (patch = model->getPatch(p->patch)))
      if (abs(p->ldim)+1 == patch->getNoParamDim())
        patch->getBoundaryNodes(abs(p->lindx),glbNodes);

  return force.initNodeMap(glbNodes,model->getNoSpaceDim());
}
