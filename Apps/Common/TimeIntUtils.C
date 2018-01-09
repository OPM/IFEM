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


TimeIntegration::Method TimeIntegration::get (const std::string& type)
{
  if (type == "be")
    return BE;
  else if (type == "bdf2")
    return BDF2;
  else if (type == "cn")
    return THETA;
  else if (type == "euler")
    return EULER;
  else if (type == "heun")
    return HEUN;
  else if (type == "heuneuler")
    return HEUNEULER;
  else if (type == "bs")
    return BOGACKISHAMPINE;
  else if (type == "fehlberg")
    return FEHLBERG;
  else if (type == "rk3")
    return RK3;
  else if (type == "rk4")
    return RK4;
  else
    return NONE;
}


int TimeIntegration::Order (Method method)
{
  switch (method) {
  case EULER:
  case BE:
    return 1;
  case HEUN:
  case BDF2:
  case THETA:
    return 2;
  case RK3:
    return 3;
  case RK4:
    return 4;
  default:
    return 0;
  }
}


int TimeIntegration::Steps (Method method)
{
  switch (method) {
  case BDF2:
    return 2;
  default:
    return 1;
  }
}
