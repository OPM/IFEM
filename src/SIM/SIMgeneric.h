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
  in a more flexible way.
*/

class SIMgeneric : public SIMoutput
{
protected:
  //! \brief Default constructor.
  explicit SIMgeneric(IntegrandBase* itg = nullptr) : SIMoutput(itg) {}

public:
  //! \brief Empty destructor.
  virtual ~SIMgeneric() {}

  //! \brief Creates a model with the default geometry (line, plane, cube).
  //! \return Pointer to the (first) spline patch of the model
  ASMbase* createDefaultModel();

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] par Parameters of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index containing the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, const double* par,
                     int deriv = 0, int patch = 1) const;

  //! \brief Evaluates the mapping of the geometry at the given point.
  //! \param[in] xi Dimensionless parameters in range [0,1] of the point
  //! \param[out] X The Cartesian coordinates of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[in] patch 1-based patch index containing the evaluation point
  //! \param[in] global If \e true, return global number, otherwise patch-local
  //! \return Patch-local or global node number of node that matches the point
  int evalPoint(const double* xi, Vec3& X, double* param = nullptr,
                int patch = 1, bool global = false) const;

  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameters of the point in the knot-span domain
  //! \param[in] patch 1-based patch index containing the point
  //! \param[in] global If \e true, return global number, otherwise patch-local
  //! \return Patch-local or global number of the element containing the point
  int findElementContaining(const double* param,
                            int patch = 1, bool global = false) const;
};

#endif
