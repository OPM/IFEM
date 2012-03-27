// $Id$
//==============================================================================
//!
//! \file NLoptions.C
//!
//! \date Mar 22 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for encapsulation of nonlinear formulation options.
//!
//==============================================================================

#include "NLoptions.h"
#include "NonlinearElasticityFbar.h"
#include "NonlinearElasticityULMixed.h"
#include "NonlinearElasticityULMX.h"
#include "NonlinearElasticity.h"
#include "SIMLinEl2D.h"
#include "ASMmxBase.h"
#include "Utilities.h"
#include "tinyxml.h"


void NLoptions::parse (const TiXmlElement* elem)
{
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"planestrain"))
      SIMLinEl2D::planeStrain = ndim == 2;
    else if (!strcasecmp(child->Value(),"axisymmetric"))
      SIMLinEl2D::axiSymmetry = ndim == 2;
    else if (!strcasecmp(child->Value(),"totallagrange"))
      form = SIM::TOTAL_LAGRANGE;
    else if (!strcasecmp(child->Value(),"updatedlagrange")) {
      if (form < SIM::UPDATED_LAGRANGE)
	form = SIM::UPDATED_LAGRANGE;
    }
    else if (!strcasecmp(child->Value(),"mixed")) {
      if (child->FirstChild())
	pOrd = atoi(child->FirstChild()->Value());
      std::string type;
      utl::getAttribute(child,"type",type);
      if (type == "Qp/Pp-1")
	form = SIM::MIXED_QnPn1;
      else if (type == "Qp/Qp-1") {
	form = SIM::MIXED_QnQn1;
	ASMmxBase::useCpminus1 = pOrd == 1;
      }
      else if (type == "Fbar")
	form = SIM::FBAR;
      else if (type == "Tensor")
	form = SIM::NONLINEAR;
    }
}


Elasticity* NLoptions::getIntegrand () const
{
  switch (form)
    {
    case SIM::FBAR:
      // F-bar formulation
      return new NonlinearElasticityFbar(ndim,SIMLinEl2D::axiSymmetry,pOrd);

    case SIM::MIXED_QnQn1:
      // Continuous volumetric change and pressure fields
      return new NonlinearElasticityULMixed(ndim,SIMLinEl2D::axiSymmetry);

    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      return new NonlinearElasticityULMX(ndim,SIMLinEl2D::axiSymmetry,pOrd);

    case SIM::UPDATED_LAGRANGE:
      return new NonlinearElasticityUL(ndim,SIMLinEl2D::axiSymmetry);

    case SIM::TOTAL_LAGRANGE:
      return new NonlinearElasticityTL(ndim,SIMLinEl2D::axiSymmetry);

    case SIM::NONLINEAR: // Old tensor-based TL-formulation
      return new NonlinearElasticity(ndim);

    default:
      std::cerr <<" *** NLoptions::getIntegrand: Unknown problem formulation "
		<< form << std::endl;
    }

  return NULL;
}
