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

#include "TimeIntUtils.h"

namespace TimeIntegration {

int Order(Method method)
{
  if (method == EULER || method == BE)
    return 1;

  if (method == HEUN || method == BDF2 || method == THETA)
    return 2;

  if (method == RK3)
    return 3;

  if (method == RK4)
    return 4;

  return 0;
}


int Steps(Method method)
{
  if (method == BDF2)
    return 2;

  return 1;
}

}
