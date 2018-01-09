//==============================================================================
//!
//! \file TimeIntUtils.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various helpers for time integration.
//!
//==============================================================================

#ifndef _TIME_INT_UTILS_H
#define _TIME_INT_UTILS_H

#include "MatVec.h"
#include <string>


namespace TimeIntegration
{
  //! \brief Enum defining various solution methods.
  //!
  //! The BE and HEUNEULER methods are used as starting markers for
  //! implicit methods and embedded methods, respecively.
  //! Therefore keep them first in these groups.
  enum Method
  {
    NONE            = 0, //!< No time integration

    // Explicit methods
    EULER           = 1, //!< Forward Euler, explicit
    HEUN            = 2, //!< Heun-Euler, explicit
    RK3             = 3, //!< Kutta's third order method, explicit
    RK4             = 4, //!< Kutta's fourth order method, explicit

    // Implicit methods
    BE              = 5, //!< Backward Euler, implicit
    BDF2            = 6, //!< Second order backward differencing, implicit

    // Embedded, explicit methods
    HEUNEULER       = 7, //!< Heun-Euler embedded order 1(2)
    BOGACKISHAMPINE = 8, //!< Bogacki-Shampine order 2(3)
    FEHLBERG        = 9, //!< Runge-Kutta-Fehlberg order 4(5)

    THETA           = 10 //!< Theta rule (includes EULER, BE and Crank-Nicolson)
  };

  //! \brief Struct holding a Runge-Kutta tableaux.
  struct RKTableaux {
    int   order; //!< Order of scheme
    Matrix    A; //!< Coefficient matrix
    RealArray b; //!< Stage weights
    RealArray c; //!< Stage levels
  };

  //! \brief Maps a text string into a Method enum value.
  Method get(const std::string& type);
  //! \brief Returns the temporal order of the given method.
  int Order(Method method);
  //! \brief Returns the number of steps the given method.
  int Steps(Method method);
}

#endif
