// $Id$
//==============================================================================
//!
//! \file NewmarkMats.h
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#ifndef _NEWMARK_MATS_H
#define _NEWMARK_MATS_H

#include "ElmMats.h"


/*!
  \brief Class representing the element matrices for a dynamic FEM problem
  based on Newmark time integration.
*/

class NewmarkMats : public ElmMats
{
public:
  //! \brief The constructor initializes the time integration parameters.
  NewmarkMats(double a1 = 0.0, double a2 = 0.0, double b = 0.0, double c = 0.0);
  //! \brief Empty destructor.
  virtual ~NewmarkMats() {}

  //! \brief Updates the time step size and the \a isPredictor flag.
  void setStepSize(double dt, int iter) { h = dt; isPredictor = iter == 0; }

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;

protected:
  double alpha1; //!< Mass-proportional damping
  double alpha2; //!< Stiffness-proportional damping
  double beta;   //!< Time integration parameter
  double gamma;  //!< Time integration parameter
  double h;      //!< Time step size

  bool isPredictor; //!< If \e true, we are in the predictor step
};

#endif
