// $Id$
//==============================================================================
//!
//! \file SIMgeneric.h
//!
//! \date Aug 28 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generic SIM class with some added functionalities.
//!
//==============================================================================

#ifndef _SIM_GENERIC_H_
#define _SIM_GENERIC_H_

#include "SIMoutput.h"


/*!
  \brief Generic SIM class with some added functionalities.
  \details This class extends the SIMbase class with some added functionalities
  of generic character, which can be used to access the FE data and structures
  on a more flexible way.
*/

class SIMgeneric : public SIMoutput
{
protected:
  //! \brief Default constructor.
  SIMgeneric(IntegrandBase* itg = NULL) : SIMoutput(itg) {}

public:
  //! \brief Emptry destructor.
  virtual ~SIMgeneric() {}

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] par Parameters of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index contining the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, const double* par,
                     int deriv = 0, int patch = 1) const;
};

#endif
