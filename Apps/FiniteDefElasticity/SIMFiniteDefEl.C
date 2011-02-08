// $Id: SIMFiniteDefEl.C,v 1.3 2010-12-29 18:51:53 kmo Exp $
//==============================================================================
//!
//! \file SIMFiniteDefEl.C
//!
//! \date Dec 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#include "SIMFiniteDefEl.h"
#include "NonlinearElasticityULMixed.h"
#include "NonlinearElasticityULMX.h"
#include "NeoHookeElasticity.h"


SIMFiniteDefEl2D::SIMFiniteDefEl2D (int form, bool planeStress,
				    const std::vector<int>& options)
  : SIMLinEl2D(SIM::NONE)
{
  int mVER = options.size() > 0 ? options[0] : -1; // material version
  int pOrd = options.size() > 1 ? options[1] :  0; // pressure field order

  switch (form)
    {
    case SIM::MIXED_QnQn1:
      nf[1] = 2; // continuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMixed(2,mVER);
      break;
    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMX(2,mVER,pOrd);
      break;
    case SIM::UPDATED_LAGRANGE:
      myProblem = new NonlinearElasticityUL(2,planeStress,mVER);
      break;
    case SIM::TOTAL_LAGRANGE:
      myProblem = new NonlinearElasticityTL(2,planeStress);
      break;
    case SIM::NONLINEAR: // Old tensor-based TL-formulation
      myProblem = new NonlinearElasticity(2,planeStress);
      break;
    default:
      std::cerr <<" *** SIMFiniteDefEl2D: Unknown problem formulation "
		<< form << std::endl;
    }
}


SIMFiniteDefEl3D::SIMFiniteDefEl3D (bool checkRHS, int form,
				    const std::vector<int>& options)
  : SIMLinEl3D(checkRHS,SIM::NONE)
{
  int mVER = options.size() > 0 ? options[0] : -1; // material version
  int pOrd = options.size() > 1 ? options[1] :  0; // pressure field order

  switch (form)
    {
    case SIM::MIXED_QnQn1:
      nf[1] = 2; // continuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMixed(3,mVER);
      break;
    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      myProblem = new NonlinearElasticityULMX(3,mVER,pOrd);
      break;
    case SIM::UPDATED_LAGRANGE:
      myProblem = new NonlinearElasticityUL(3,false,mVER);
      break;
    case SIM::TOTAL_LAGRANGE:
      myProblem = new NonlinearElasticityTL();
      break;
    case SIM::NONLINEAR: // Old tensor-based TL-formulation
      myProblem = new NonlinearElasticity();
      break;
    case SIM::NEOHOOKE: // NeoHookean TL-formulation (not working)
      myProblem = new NeoHookeElasticity();
      break;
    case SIM::NEOHOOKE_IV: // NeoHookean isochoric/volumetric TL-formulation
      myProblem = new NeoHookeElasticityIV();
      break;
    default:
      std::cerr <<" *** SIMFiniteDefEl3D: Unknown problem formulation "
		<< form << std::endl;
    }
}
