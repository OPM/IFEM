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
  //! \param[in] a1 Mass-proportional damping coefficient
  //! \param[in] a2 Stiffness-proportional damping coefficient
  //! \param[in] b Time integration parameter
  //! \param[in] c Time integration parameter
  //! \param[in] generalizedAlpha If \e true, interpret \a b and \a c as the
  //! generalized-alpha parameters &alpha;<sub>m</sub> and &alpha;<sub>f</sub>,
  //! respectively, otherwise as &beta; and &gamma;, respectively
  NewmarkMats(double a1, double a2, double b = 0.0, double c = 0.0,
              bool generalizedAlpha = false);
  //! \brief Empty destructor.
  virtual ~NewmarkMats() {}

  //! \brief Updates the time step size and the \ref isPredictor flag.
  //! \param[in] dt New time step size
  //! \param[in] it Newton-Raphson iteration counter
  virtual void setStepSize(double dt, int it) { h = dt; isPredictor = it == 0; }

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;

protected:
  bool  isPredictor; //!< If \e true, we are in the predictor step
  double h;          //!< Time step size

  double alpha1; //!< Mass-proportional damping coefficient
  double alpha2; //!< Stiffness-proportional damping coefficient
  double beta;   //!< Newmark time integration parameter, &beta;
  double gamma;  //!< Newmark time integration parameter, &gamma;

private:
  double alpha_m; //!< Generalized-alpha parameter, &alpha;<sub>m</sub>
  double alpha_f; //!< Generalized-alpha parameter, &alpha;<sub>f</sub>
  bool   slvDisp; //!< If \e true, solve for displacement increments
};

#endif
