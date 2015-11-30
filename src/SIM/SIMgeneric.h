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
  //! \brief Empty destructor.
  virtual ~SIMgeneric() {}

  //! \brief Creates a model with the default geometry (line, plane, cube).
  void createDefaultModel();

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] par Parameters of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index contining the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, const double* par,
                     int deriv = 0, int patch = 1) const;

  //! \brief Evaluates the mapping of the geometry at the given point.
  //! \param[in] xi Dimensionless parameters in range [0,1] of the point
  //! \param[out] X The Cartesian coordinates of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[in] patch The patch to evaluate
  //! \return 0 if the evaluation went good
  int evalPoint(const double* xi, Vec3& X, 
                double* param=NULL, int patch = 1) const;
};

#endif
