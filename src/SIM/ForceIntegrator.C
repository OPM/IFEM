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
#include "Profiler.h"


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
  PROFILE1("SIM::getBoundaryForce");

  ForceBase* forceInt = model->getBoundaryForceIntegrand(X0);
  if (!forceInt)
  {
    std::cerr <<" *** SIM::getBoundaryForce: No force integrand."<< std::endl;
    return Vector();
  }

  forceInt->initBuffer(model->getNoElms());

  Vector force;
  if (integrate(solution,model,code,time,forceInt))
    forceInt->assemble(force); // Assemble into the global force resultant
  else
    std::cerr <<" *** SIM::getBoundaryForce: Failed to evaluate boundary force"
              << std::endl;

  delete forceInt;
  return force;
}


bool SIM::getNodalForces (const Vectors& solution, SIMbase* model, int code,
                          const TimeDomain& time, GlbForceVec& force)
{
  PROFILE1("SIM::getNodalForces");

  force.initialize();

  bool ok = false;
  ForceBase* forceInt = model->getNodalForceIntegrand();
  if (!forceInt)
    std::cerr <<" *** SIM::getNodalForces: No force integrand."<< std::endl;
  else if (!integrate(solution,model,code,time,forceInt,&force))
    std::cerr <<" *** SIM::getNodalForces: Failed to evaluate nodal forces"
              << std::endl;
  else
    ok = force.finalize();

  delete forceInt;
  return ok;
}


bool SIM::initBoundaryNodeMap (SIMbase* model, int code, GlbForceVec& force)
{
  ASMbase* patch;
  IntVec glbNodes;
  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop(); p++)
    if (abs(p->pindx) == code && (patch = model->getPatch(p->patch)))
      if (abs(p->ldim)+1 == patch->getNoParamDim())
        patch->getBoundaryNodes(abs(p->lindx)%10,glbNodes);

  return force.initNodeMap(glbNodes,model->getNoSpaceDim());
}


bool SIM::integrate(const Vectors& solution, SIMbase* model, int code,
                    const TimeDomain& time, ForceBase* forceInt,
                    GlbForceVec* force)
{
  // Integrate forces for given boundary segment
  bool ok = true;
  size_t prevPatch = 0;
  GlobalIntegral dummy;
  GlobalIntegral& frc = force ? *force : dummy;

  if (code == 0) { // special case - volume integral over the entire model
    for (int p = 1; p <= model->getNoPatches() && ok; ++p) {
      ASMbase* patch = model->getPatch(p);
      ok = model->extractPatchSolution(solution,p-1);
      model->setPatchMaterial(p);
      ok &= patch->integrate(*forceInt,frc,time);
    }
    return ok;
  }

  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop() && ok; p++)
    if (abs(p->pindx) == code)
    {
      ASMbase* patch = model->getPatch(p->patch);
      if (!patch)
        ok = false;
      else if ( forceInt->hasInteriorTerms() ||
               (forceInt->hasBoundaryTerms() &&
                abs(p->ldim)+1 == patch->getNoParamDim()))
      {
        if (p->patch != prevPatch) {
          ok = model->extractPatchSolution(solution,p->patch-1);
          model->setPatchMaterial(p->patch);
        }
        if (forceInt->hasInteriorTerms())
          ok &= patch->integrate(*forceInt,frc,time);
        else
          ok &= patch->integrate(*forceInt,abs(p->lindx),frc,time);
        prevPatch = p->patch;
      }
    }

  return ok;
}
