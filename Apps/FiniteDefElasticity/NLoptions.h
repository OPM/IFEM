// $Id$
//==============================================================================
//!
//! \file NLoptions.h
//!
//! \date Mar 22 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for encapsulation of nonlinear formulation options.
//!
//==============================================================================

#ifndef _NL_OPTIONS_H
#define _NL_OPTIONS_H

#include "SIMenums.h"

class TiXmlElement;
class Elasticity;


namespace SIM
{
  //! \brief Enum defining various finite deformation formulations.
  enum NlFormulation
  {
    TOTAL_LAGRANGE   = 3,
    UPDATED_LAGRANGE = 4,
    MIXED_QnPn1      = 5,
    MIXED_QnQn1      = 6,
    FBAR             = 7
  };
};


/*!
  \brief Class for encapsulation of nonlinear formulation options.
*/

class NLoptions
{
public:
  //! \brief Default constructor.
  NLoptions(int n = 0) : ndim(n), form(SIM::TOTAL_LAGRANGE), pOrd(0) {}

  //! \brief Parses input options from an XML-tag.
  void parse(const TiXmlElement* elem);

  //! \brief Returns a dynamically allocated integrand for given formulation.
  Elasticity* getIntegrand() const;

public:
  int ndim; //!< Number of parametric dimensions (2 or 3)
  int form; //!< Nonlinear formulation option
  int pOrd; //!< Order of internal pressure field
};

#endif
