//==============================================================================
//!
//! \file TimeIntUtils.C
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various helpers for time integration.
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
  else if (type == "ab1")
    return AB1;
  else if (type == "ab2")
    return AB2;
  else if (type == "ab3")
    return AB3;
  else if (type == "ab4")
    return AB4;
  else if (type == "ab5")
    return AB5;
  else if (type == "am1")
    return AM1;
  else if (type == "am2")
    return AM2;
  else if (type == "am3")
    return AM3;
  else if (type == "am4")
    return AM4;
  else
    return NONE;
}


int TimeIntegration::Order (Method method)
{
  switch (method) {
  case AB1:
  case AM1:
  case BE:
  case EULER:
    return 1;
  case AB2:
  case AM2:
  case BDF2:
  case HEUN:
  case THETA:
    return 2;
  case AB3:
  case AM3:
  case RK3:
    return 3;
  case AB4:
  case AM4:
  case RK4:
    return 4;
  case AB5:
    return 5;
  default:
    return 0;
  }
}


int TimeIntegration::Steps (Method method)
{
  switch (method) {
  case AB5:
    return 5;
  case AB4:
  case AM4:
    return 4;
  case AB3:
  case AM3:
    return 3;
  case AB2:
  case AM2:
  case BDF2:
    return 2;
  default:
    return 1;
  }
}
