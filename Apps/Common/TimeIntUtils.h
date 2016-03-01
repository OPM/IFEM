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


namespace TimeIntegration //! Time integration scope
{
  //! \brief Enum defining various solution methods.
  // The BE and HEUNEULER methods are used as starting markers for implicit
  // methods and embedded methods. Keep them first in these groups
  enum Method
  {
    NONE             = 0,  //!< No time integration

    // Explicit methods
    EULER            = 1,  //!< Forward Euler, explicit
    HEUN             = 2,  //!< Heun-Euler, explicit
    RK3              = 3,  //!< Kutta's third order method, explicit
    RK4              = 4,  //!< Kutta's fourth order method, explicit

    // Implicit methods
    BE               = 5,  //!< Backward Euler, implicit
    BDF2             = 6,  //!< Second order backward differencing, implicit

    // Embedded, explicit methods
    HEUNEULER        = 7,  //!< Heun-Euler embedded order 1(2)
    BOGACKISHAMPINE  = 8,  //!< Bogacki-Shampine order 2(3)
    FEHLBERG         = 9,  //!< Runge-Kutta-Fehlberg order 4(5)

    THETA            = 10  //!< Theta rule (family includes EULER, BE and Crank-Nicolson)
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
