// $Id$
//==============================================================================
//!
//! \file ForceIntegrator.C
//!
//! \date May 10 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for integration of boundary forces.
//!
//==============================================================================

#include "ForceIntegrator.h"
#include "SIMbase.h"
#include "ASMbase.h"
#include "IntegrandBase.h"
#include "GlobalIntegral.h"
#ifdef PARALLEL_PETSC
#include "petscversion.h"
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
#include "petscpcmg.h"
#else
#include "petscmg.h"
#endif
#include <mpi.h>
#endif


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

  const std::vector<ASMbase*>& feModel = model->getFEModel();

  // Integrate forces for given boundary segment
  bool ok = true;
  size_t prevPatch = 0;
  GlobalIntegral dummy;
  PropertyVec::const_iterator p;
  for (p = model->begin_prop(); p != model->end_prop() && ok; p++)
    if (abs(p->pindx) == code)
    {
      size_t j = p->patch;
      if (j < 1 || j > feModel.size())
        ok = false;
      else if (abs(p->ldim)+1 == feModel[j-1]->getNoSpaceDim())
      {
        if (j != prevPatch)
          ok = model->extractPatchSolution(solution,j-1);
        ok &= feModel[j-1]->integrate(*forceInt,abs(p->lindx),dummy,time);
        prevPatch = j;
      }
    }

  // Assemble the element force contributions into the global force resultant
  Vector force;
  forceInt->assemble(force);
  delete forceInt;

#ifdef PARALLEL_PETSC
  size_t nComp = force.size();
  double* tmp = new double[nComp];
  MPI_Allreduce(force.ptr(),tmp,nComp,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  memcpy(force.ptr(),tmp,nComp*sizeof(double));
  delete[] tmp;
#endif

  if (ok) return force;

  std::cerr <<" *** SIM::getBoundaryForce: Failed to evaluate boundary force"
	    << std::endl;
  return Vector();
}
