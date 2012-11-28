//==============================================================================
//!
//! \file TimeIntUtils.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various helpers for time integration
//!
//==============================================================================

#ifndef _TIME_INT_UTILS_H
#define _TIME_INT_UTILS_H

#include "DenseMatrix.h"
#include <vector>


namespace TimeIntegration //!< Time integration scope
{
  //! \brief Enum defining various solution methods.
  enum Method
  {
    NONE  = 0,  //!< No time integration
    EULER = 1,  //!< Forward Euler, explicit
    HEUN  = 2,  //!< Heun-Euler, explicit
    RK3   = 3,  //!< Kutta's third order method, explicit
    RK4   = 4,  //!< Kutta's fourth order method, explicit
    BE    = 5,  //!< Backward Euler, implicit
    BDF2  = 6   //!< Second order backward differencing, implicit
  };

  //! \brief Struct holding a Runge-Kutta tableaux
  typedef struct {
    int order;             //!< order of Scheme
    DenseMatrix A;         //!< Coefficient matrix
    std::vector<double> b; //!< Stage weights
    std::vector<double> c; //!< Stage levels
  } RKTableaux;

  //! \brief Returns the temporal order of the given method
  //! \param[in] method The method
  //! \returns The order of the given method
  int Order(Method method);

  //! \brief Returns the number of steps the given method
  //! \param[in] method The method
  //! \returns The number of steps in the given method
  int Steps(Method method);
}

#endif
